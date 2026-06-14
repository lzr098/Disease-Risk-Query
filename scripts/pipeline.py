"""End-to-end disease risk query pipeline."""

from __future__ import annotations

import json
import logging
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from constants import (
    AD_CORE_GENES,
    AD_GWAS_LOCI_GENES,
    AD_MENDELIAN_CORE_GENES,
    DEFAULT_DISEASE_MODE,
    DEFAULT_OUTPUT_DIR,
    GRCH38_FASTA,
    resolve_builtin_disease_key,
    resolve_disease_mode,
)
from apoe_checker import check_apoe
from clinvar_phenotype_matcher import filter_variants_by_clinvar_disease
from disease_reference import get_disease_reference
from gene_set_builder import build_disease_gene_set
from gpa_runner import run_gpa_on_filtered_vcf
from gwas_source import annotate_gwas_hits, check_gwas_lead_snps, score_gwas_support
from hpo_mapper import resolve_disease_query
from literature_source import annotate_literature_support, score_literature_support
from liftover import detect_genome_build, liftover_vcf, validate_ref_alleles
from report import generate_report
from risk_scorer import calculate_contribution_score, calculate_total_score, sex_age_bonus
from tier_filters import denoise_tier3, remove_duplicate_variants
from vcf_filter import check_vcf_completeness, filter_vcf_by_gene_set

logger = logging.getLogger(__name__)


@dataclass
class PipelineConfig:
    input_vcf: Path
    disease_query: str
    hpo_id: Optional[str] = None
    sex: str = "unknown"
    age: Optional[int] = None
    family_history: bool = False
    tissue: Optional[str] = None
    output_dir: Path = DEFAULT_OUTPUT_DIR
    max_genes: int = 200
    offline: bool = False
    spliceai: bool = True
    two_phase: bool = False
    skip_liftover: bool = False
    chain_file: Optional[Path] = None
    literature_genes: Optional[list[str]] = None
    literature_variants: Optional[list[dict]] = None
    refresh_cache: bool = False
    disease_mode: str = DEFAULT_DISEASE_MODE


def run_disease_risk_pipeline(config: PipelineConfig) -> dict:
    """Execute full Step 1-6 pipeline and return results + report path."""
    output_dir = config.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    disease_mode = resolve_disease_mode(config.disease_query, config.disease_mode)
    logger.info("Disease mode resolved to: %s", disease_mode)

    run_id = f"drq_{config.disease_query.replace(' ', '_').replace('/', '_')[:30]}"
    work_dir = output_dir / run_id
    work_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Genome build normalization
    logger.info("Step 1: Detecting genome build")
    build = detect_genome_build(config.input_vcf)
    logger.info("Detected build: %s", build)

    normalized_vcf = work_dir / "normalized.vcf.gz"
    if build in {"GRCh37", "hg19"} and not config.skip_liftover:
        lo_result = liftover_vcf(
            config.input_vcf,
            normalized_vcf,
            source_build=build,
            chain_file=config.chain_file,
        )
        normalized_vcf = Path(lo_result["output_path"])
        logger.info("Liftover result: %s", lo_result)
        build = "GRCh38"
    else:
        shutil.copy2(config.input_vcf, normalized_vcf)
        idx = config.input_vcf.with_suffix(config.input_vcf.suffix + ".csi")
        if idx.exists():
            shutil.copy2(idx, normalized_vcf.with_suffix(normalized_vcf.suffix + ".csi"))

    # REF validation
    if GRCH38_FASTA.exists() and build == "GRCh38":
        logger.info("Step 1b: Validating REF alleles against GRCh38 FASTA")
        validation = validate_ref_alleles(normalized_vcf, GRCH38_FASTA)
        logger.info("REF validation: %s", validation)
        if not validation["pass"]:
            raise RuntimeError(
                f"REF validation failed: mismatch rate {validation['mismatch_rate']:.1%}. "
                "Check input VCF integrity or genome build."
            )
    else:
        validation = {"pass": True, "note": "FASTA not available or build not GRCh38"}

    # Step 1c: VCF completeness / filtering detection
    logger.info("Step 1c: Checking VCF completeness (common variant filtering)")
    vcf_qc = check_vcf_completeness(normalized_vcf, config.disease_query)
    logger.info("VCF QC: %s", vcf_qc)

    # Step 2b: APOE genotype check (Alzheimer-specific)
    apoe_result = None
    canonical_key = resolve_builtin_disease_key(config.disease_query)
    if canonical_key == "alzheimer disease" or "alzheimer" in config.disease_query.lower():
        logger.info("Step 2b: Checking APOE genotype")
        apoe_result = check_apoe(normalized_vcf)
        logger.info("APOE result: %s", apoe_result)

    # Step 2c: Direct check of known GWAS lead SNPs
    logger.info("Step 2c: Checking known GWAS lead SNPs")
    gwas_lead_snp_results = check_gwas_lead_snps(normalized_vcf, config.disease_query)
    gwas_lead_hits = [s for s in gwas_lead_snp_results if s.get("sample_gt")]
    logger.info("GWAS lead SNPs checked: %d, hits: %d", len(gwas_lead_snp_results), len(gwas_lead_hits))

    # Step 3: Disease -> HPO
    logger.info("Step 3: Mapping disease query to HPO")
    hpo_result = resolve_disease_query(config.disease_query, config.hpo_id)
    hpo_id = hpo_result.get("hpo_id")
    hpo_name = hpo_result.get("hpo_name")
    if hpo_id is None:
        logger.warning("No HPO mapping found for '%s'; proceeding with OMIM-only gene set", config.disease_query)

    # Step 3b: Load or build disease reference cache
    logger.info("Step 3b: Loading disease reference cache")
    disease_ref = get_disease_reference(config.disease_query, refresh=config.refresh_cache)
    ref_core_genes = set(disease_ref.get("core_genes", []))
    gene_contribution_map = disease_ref.get("gene_contribution_map", {})
    snp_contribution_map = disease_ref.get("snp_contribution_map", {})
    logger.info(
        "Disease reference: %d core genes, %d literature entries, %d gene contributions, %d snp contributions",
        len(ref_core_genes), len(disease_ref.get("literature", [])),
        len(gene_contribution_map), len(snp_contribution_map),
    )

    # Step 4: Build disease gene set
    logger.info("Step 4: Building disease gene set")
    ref_gwas_loci = set(disease_ref.get("gwas_loci_genes", []))
    core_genes = ref_core_genes if ref_core_genes else (AD_CORE_GENES if "alzheimer" in config.disease_query.lower() else set())
    # Use Mendelian core for strict tier filtering; GWAS loci are informative but not relaxed
    mendelian_core = ref_core_genes if ref_core_genes else (AD_MENDELIAN_CORE_GENES if "alzheimer" in config.disease_query.lower() else set())
    # Core + GWAS loci from reference serve as literature-backed extra genes
    extra_genes = sorted(set(config.literature_genes or []) | core_genes | ref_gwas_loci)
    # OMIM titles/symbols are English-only; use the canonical English disease
    # key for OMIM keyword search so Chinese queries still find OMIM entries.
    omim_search_term = canonical_key or config.disease_query
    gene_set_result = build_disease_gene_set(
        disease_name=config.disease_query,
        hpo_genes=hpo_result.get("genes") if hpo_id else [],
        omim_keywords=[omim_search_term],
        extra_genes=extra_genes,
        max_genes=config.max_genes,
        core_genes=core_genes,
    )
    gene_set_result["disease_reference"] = disease_ref

    if not gene_set_result["merged_genes"]:
        # No genes found; produce empty report
        report_path = work_dir / "report.md"
        empty_gpa = {
            "tier1_variants": [],
            "tier2_variants": [],
            "tier3_variants": [],
            "multi_hit": [],
            "summary": {},
        }
        empty_score = calculate_total_score([], [], [], [], set(), 0)
        generate_report(
            disease_name=config.disease_query,
            hpo_id=hpo_id,
            hpo_name=hpo_name,
            sex=config.sex,
            age=config.age,
            gene_set_result=gene_set_result,
            gpa_result=empty_gpa,
            score_result=empty_score,
            output_path=report_path,
            apoe_result=apoe_result,
            vcf_qc=vcf_qc,
            disease_mode=disease_mode,
        )
        return {
            "success": True,
            "hpo_id": hpo_id,
            "hpo_name": hpo_name,
            "gene_set": gene_set_result,
            "gpa": empty_gpa,
            "score": empty_score,
            "apoe": apoe_result,
            "report_path": report_path,
            "work_dir": work_dir,
            "validation": validation,
            "note": "No disease-associated genes found in local databases.",
        }

    # Step 5A: Filter VCF to disease gene regions
    logger.info("Step 5A: Filtering VCF to %d disease genes", gene_set_result["total"])
    filtered_vcf = work_dir / "filtered_disease_genes.vcf.gz"
    filter_stats = filter_vcf_by_gene_set(
        normalized_vcf,
        filtered_vcf,
        gene_set_result["merged_genes"],
        clinvar_safetynet=True,
        disease_name=config.disease_query,
    )
    logger.info("Filter stats: %s", filter_stats)

    # Step 5B: GPA tier classification
    logger.info("Step 5B: Running GPA tier classification")
    progress_log = work_dir / "gpa_progress.jsonl"
    gpa_result = run_gpa_on_filtered_vcf(
        filtered_vcf=filtered_vcf,
        disease_name=config.disease_query,
        tissue=config.tissue,
        sex=config.sex,
        age=config.age,
        offline=config.offline,
        spliceai=config.spliceai,
        two_phase=config.two_phase,
        progress_log=progress_log,
    )

    if not gpa_result.get("success"):
        return {
            "success": False,
            "error": gpa_result.get("error", "GPA analysis failed"),
            "hpo_id": hpo_id,
            "gene_set": gene_set_result,
            "work_dir": work_dir,
        }

    # Step 5C: Post-GPA filtering
    logger.info("Step 5C: Post-GPA ClinVar phenotype matching and Tier 3 denoising")
    # Tier 1 variants should also be required to match the disease phenotype in
    # ClinVar; otherwise broad pathogenic variants in non-disease genes (e.g.
    # VWF p.Gln1311Ter) are miscalled as Tier 1.
    tier1_before = len(gpa_result.get("tier1_variants", []))
    gpa_result["tier1_variants"] = filter_variants_by_clinvar_disease(
        gpa_result.get("tier1_variants", []),
        config.disease_query,
        require_match=True,
    )
    logger.info(
        "Tier 1 ClinVar phenotype filtering: %d -> %d",
        tier1_before, len(gpa_result["tier1_variants"]),
    )
    tier2_before = len(gpa_result.get("tier2_variants", []))
    gpa_result["tier2_variants"] = filter_variants_by_clinvar_disease(
        gpa_result.get("tier2_variants", []),
        config.disease_query,
        require_match=True,
    )
    logger.info(
        "Tier 2 ClinVar phenotype filtering: %d -> %d",
        tier2_before, len(gpa_result["tier2_variants"]),
    )

    # Deduplicate and denoise Tier 3
    gpa_result["tier3_variants"] = remove_duplicate_variants(
        gpa_result.get("tier3_variants", [])
    )
    gpa_result["tier3_variants"] = denoise_tier3(
        gpa_result["tier3_variants"],
        core_genes=mendelian_core or AD_MENDELIAN_CORE_GENES,
        gwas_loci_genes=ref_gwas_loci or AD_GWAS_LOCI_GENES,
    )

    # Rebuild multi-hit from filtered tiers
    all_tiered = (
        gpa_result.get("tier1_variants", [])
        + gpa_result.get("tier2_variants", [])
        + gpa_result.get("tier3_variants", [])
    )
    from collections import Counter
    gene_counts = Counter(v.get("gene") or v.get("GENE") for v in all_tiered if (v.get("gene") or v.get("GENE")))
    gpa_result["multi_hit"] = [
        {"gene": gene, "variant_count": count}
        for gene, count in gene_counts.items()
        if count > 1
    ]

    # Step 5D: Annotate tiers with GWAS and literature evidence
    logger.info("Step 5D: Annotating variants with GWAS and literature evidence")
    annotate_gwas_hits(all_tiered, config.disease_query)
    annotate_literature_support(all_tiered, config.disease_query)
    gwas_summary = score_gwas_support(all_tiered, config.disease_query)
    lit_summary = score_literature_support(all_tiered, config.disease_query, core_genes=mendelian_core)
    logger.info("GWAS proximal variants: %d in genes %s", gwas_summary["hit_count"], gwas_summary["hit_genes"])
    logger.info("Literature-supported variants: %d in genes %s", len(lit_summary["variant_hits"]), lit_summary["gene_hits"])

    # Step 6: Composite risk scoring
    logger.info("Step 6: Calculating composite %s score", disease_mode)
    if disease_mode == "complex":
        score_result = calculate_contribution_score(
            tier1=gpa_result.get("tier1_variants", []),
            tier2=gpa_result.get("tier2_variants", []),
            tier3=gpa_result.get("tier3_variants", []),
            core_genes=ref_core_genes,
            gwas_lead_snps=gwas_lead_snp_results,
            literature_variants=config.literature_variants,
            literature_genes=set(config.literature_genes or []),
            vcf_qc=vcf_qc,
            gene_contribution_map=gene_contribution_map,
            snp_contribution_map=snp_contribution_map,
        )
    else:
        sab = sex_age_bonus(config.disease_query, config.sex, config.age, gene_set_result["merged_genes"])
        multi_hit_genes = [
            (m["gene"] if isinstance(m, dict) else m)
            for m in gpa_result.get("multi_hit", [])
        ]
        score_result = calculate_total_score(
            tier1=gpa_result.get("tier1_variants", []),
            tier2=gpa_result.get("tier2_variants", []),
            multi_hit_genes=multi_hit_genes,
            literature_variants=config.literature_variants,
            literature_genes=set(config.literature_genes or []),
            sex_age_bonus=sab,
            gene_contribution_map=gene_contribution_map,
        )

    # Apply literature/GWAS flags to variants for report
    for v in gpa_result.get("tier1_variants", []) + gpa_result.get("tier2_variants", []):
        flags = []
        if v.get("literature_support"):
            flags.append("LITERATURE_SUPPORTED")
        if v.get("gwas_proximal"):
            flags.append(f"GWAS_PROXIMAL:{v.get('nearest_gwas_snp')}")
        v["_drq_flags"] = flags

    # Generate report
    report_path = work_dir / "report.md"
    generate_report(
        disease_name=config.disease_query,
        hpo_id=hpo_id,
        hpo_name=hpo_name,
        sex=config.sex,
        age=config.age,
        gene_set_result=gene_set_result,
        gpa_result=gpa_result,
        score_result=score_result,
        apoe_result=apoe_result,
        gwas_summary=gwas_summary,
        literature_summary=lit_summary,
        disease_reference=disease_ref,
        gwas_lead_snps=gwas_lead_snp_results,
        vcf_qc=vcf_qc,
        output_path=report_path,
        disease_mode=disease_mode,
    )

    # Save structured JSON
    json_path = work_dir / "result.json"
    result = {
        "success": True,
        "hpo_id": hpo_id,
        "hpo_name": hpo_name,
        "input_vcf": str(config.input_vcf),
        "normalized_vcf": str(normalized_vcf),
        "filtered_vcf": str(filtered_vcf),
        "genome_build": build,
        "validation": validation,
        "apoe": apoe_result,
        "disease_mode": disease_mode,
        "gene_set": gene_set_result,
        "filter_stats": filter_stats,
        "gpa": gpa_result,
        "score": score_result,
        "gwas_lead_snps": gwas_lead_snp_results,
        "gwas_summary": gwas_summary,
        "literature_summary": lit_summary,
        "disease_reference": disease_ref,
        "vcf_qc": vcf_qc,
        "report_path": str(report_path),
        "work_dir": str(work_dir),
    }
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, ensure_ascii=False, default=str)

    return result
