#!/usr/bin/env python3
"""CLI entry point for gpa-disease-risk-query skill."""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
from pathlib import Path
from typing import Optional

from pipeline import PipelineConfig, run_disease_risk_pipeline


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="gpa-disease-risk-query",
        description="Query disease-specific germline genetic risk from a WGS/WES VCF.",
    )
    parser.add_argument("--vcf", "--input", type=Path, help="Input VCF/VCF.gz path")
    parser.add_argument("--disease", help="Disease name or natural-language query")
    parser.add_argument("--hpo-id", help="Explicit HPO ID (e.g. HP:0000726)")
    parser.add_argument("--sex", default="unknown", choices=["male", "female", "unknown"], help="Sample sex")
    parser.add_argument("--age", type=int, help="Sample age")
    parser.add_argument("--family-history", action="store_true", help="First-degree relative affected")
    parser.add_argument("--tissue", help="Override GPA tissue profile (general/neurological/cardiovascular/renal/hepatic/hematopoietic)")
    parser.add_argument("--output-dir", type=Path, default=Path.cwd(), help="Output directory")
    parser.add_argument("--max-genes", type=int, default=200, help="Maximum disease genes to include")
    parser.add_argument("--offline", action="store_true", help="Run GPA in offline mode")
    parser.add_argument("--no-spliceai", action="store_true", help="Disable SpliceAI")
    parser.add_argument("--two-phase", action="store_true", help="Enable two-phase GPA pipeline")
    parser.add_argument("--skip-liftover", action="store_true", help="Skip GRCh37->GRCh38 liftover even if detected")
    parser.add_argument("--chain-file", type=Path, help="Custom liftOver chain file")
    parser.add_argument("--literature-genes", help="Comma-separated literature-derived genes")
    parser.add_argument("--literature-variants", type=Path, help="JSON file with literature variants [{CHROM,POS,REF,ALT,...}]")
    parser.add_argument("--refresh-cache", action="store_true", help="Force rebuild of disease reference cache")
    parser.add_argument("--no-gwas", action="store_true", help="Disable GWAS evidence annotation")
    parser.add_argument("--no-literature", action="store_true", help="Disable literature evidence annotation")
    parser.add_argument(
        "--disease-mode",
        choices=["mendelian", "complex", "auto"],
        default="auto",
        help=(
            "Scoring/reporting mode: mendelian uses Tier 1/2/3 risk levels; "
            "complex uses genetic contribution scoring for polygenic traits; "
            "auto picks based on built-in disease template or keyword heuristics."
        ),
    )
    parser.add_argument("--preflight", action="store_true", help="Run dependency checks and exit")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    return parser


def _parse_literature_variants(path: Optional[Path]) -> Optional[list[dict]]:
    if path is None:
        return None
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    if isinstance(data, list):
        return data
    if isinstance(data, dict) and "variants" in data:
        return data["variants"]
    raise ValueError("literature-variants JSON must be a list or {variants: [...]}")


def preflight_check() -> dict:
    """Check required local resources."""
    from constants import (
        CLINVAR_PATHOGENIC_BED,
        CLINVAR_VCF,
        GENCODE_GENE_LOCI_BED,
        GENCODE_GTF,
        GRCH38_FASTA,
        HGNC_LOOKUP,
        HPO_GENES_TO_PHENOTYPE,
        OMIM_DB,
    )

    items = []
    for name, path in [
        ("GRCh38 FASTA", GRCH38_FASTA),
        ("ClinVar VCF", CLINVAR_VCF),
        ("OMIM SQLite", OMIM_DB),
        ("HPO genes->phenotype", HPO_GENES_TO_PHENOTYPE),
        ("HGNC lookup", HGNC_LOOKUP),
        ("GENCODE GTF", GENCODE_GTF),
        ("dgra-prefilter gene BED", GENCODE_GENE_LOCI_BED),
        ("dgra-prefilter ClinVar BED", CLINVAR_PATHOGENIC_BED),
    ]:
        exists = path.exists()
        items.append({
            "name": name,
            "path": str(path),
            "exists": exists,
            "status": "PASS" if exists else "FAIL",
        })

    # Check commands
    import shutil
    for cmd in ["bcftools", "samtools"]:
        found = shutil.which(cmd) is not None
        items.append({"name": cmd, "exists": found, "status": "PASS" if found else "FAIL"})

    # Check Python libraries
    for lib in ["vcfpy", "pyliftover"]:
        try:
            __import__(lib)
            items.append({"name": lib, "exists": True, "status": "PASS"})
        except Exception as exc:
            items.append({"name": lib, "exists": False, "status": f"FAIL: {exc}"})

    # Check VEP Docker image
    try:
        docker_ok = subprocess.run(
            ["docker", "image", "inspect", "ensemblorg/ensembl-vep:latest"],
            capture_output=True, check=False,
        ).returncode == 0
        items.append({"name": "VEP Docker image", "exists": docker_ok, "status": "PASS" if docker_ok else "FAIL"})
    except Exception as exc:
        items.append({"name": "VEP Docker image", "exists": False, "status": f"FAIL: {exc}"})

    # Check VEP cache (optional but strongly recommended)
    from constants import VEP_CACHE
    cache_ok = (VEP_CACHE / "homo_sapiens").exists()
    items.append({"name": "VEP cache", "exists": cache_ok, "status": "PASS" if cache_ok else "WARN"})

    # Check liftover chain availability (do not auto-download during preflight)
    from liftover import LIFTOVER_DIR
    expected_chain = LIFTOVER_DIR / "hg19ToHg38.over.chain.gz"
    chain_exists = expected_chain.exists()
    items.append({
        "name": "liftover chain (GRCh37->GRCh38)",
        "exists": chain_exists,
        "status": "PASS" if chain_exists else "FAIL (will auto-download on first liftover)",
    })

    return {
        "overall_ready": all(
            i["exists"] for i in items
            if i["name"] not in {"HGNC lookup", "VEP cache"}
        ),
        "items": items,
    }


def main(argv: Optional[list[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        stream=sys.stderr,
    )

    if args.preflight:
        report = preflight_check()
        print(json.dumps(report, indent=2, ensure_ascii=False))
        return 0 if report["overall_ready"] else 1

    if not args.vcf or not args.disease:
        parser.error("--vcf and --disease are required (unless using --preflight)")

    lit_variants = _parse_literature_variants(args.literature_variants)
    lit_genes = [g.strip() for g in args.literature_genes.split(",") if g.strip()] if args.literature_genes else None

    config = PipelineConfig(
        input_vcf=args.vcf,
        disease_query=args.disease,
        hpo_id=args.hpo_id,
        sex=args.sex,
        age=args.age,
        family_history=args.family_history,
        tissue=args.tissue,
        output_dir=args.output_dir,
        max_genes=args.max_genes,
        offline=args.offline,
        spliceai=not args.no_spliceai,
        two_phase=args.two_phase,
        skip_liftover=args.skip_liftover,
        chain_file=args.chain_file,
        literature_genes=lit_genes,
        literature_variants=lit_variants,
        refresh_cache=args.refresh_cache,
        disease_mode=args.disease_mode,
        gwas_enabled=not args.no_gwas,
        literature_enabled=not args.no_literature,
    )

    try:
        result = run_disease_risk_pipeline(config)
    except Exception as exc:
        logging.exception("Pipeline failed")
        print(json.dumps({"success": False, "error": str(exc)}, indent=2, ensure_ascii=False))
        return 1

    if not result.get("success"):
        print(json.dumps(result, indent=2, ensure_ascii=False, default=str))
        return 1

    print(json.dumps({
        "success": True,
        "hpo_id": result.get("hpo_id"),
        "hpo_name": result.get("hpo_name"),
        "gene_count": result["gene_set"]["total"],
        "gene_set_truncated": result["gene_set"].get("truncated", False),
        "tier1": len(result["gpa"].get("tier1_variants", [])),
        "tier2": len(result["gpa"].get("tier2_variants", [])),
        "tier3": len(result["gpa"].get("tier3_variants", [])),
        "total_score": result["score"]["total_score"],
        "risk_level": result["score"].get("risk_level"),
        "contribution_level": result["score"].get("contribution_level"),
        "disease_mode": result.get("disease_mode", "mendelian"),
        "apoe": result.get("apoe"),
        "vcf_common_variants_filtered": result.get("vcf_qc", {}).get("common_variants_filtered", False),
        "vcf_presence_rate": result.get("vcf_qc", {}).get("presence_rate"),
        "gwas_hits": result.get("gwas_summary", {}).get("hit_count", 0),
        "gwas_genes": result.get("gwas_summary", {}).get("hit_genes", []),
        "gwas_lead_snp_hits": len([s for s in result.get("gwas_lead_snps", []) if s.get("sample_gt")]),
        "literature_hits": len(result.get("literature_summary", {}).get("variant_hits", [])),
        "literature_genes": result.get("literature_summary", {}).get("gene_hits", []),
        "report_path": result["report_path"],
        "work_dir": result["work_dir"],
    }, indent=2, ensure_ascii=False, default=str))

    return 0


if __name__ == "__main__":
    sys.exit(main())
