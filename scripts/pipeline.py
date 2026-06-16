"""End-to-end disease risk query pipeline (DiseaseProfile-driven v0.11)."""

from __future__ import annotations

import json
import logging
import shutil
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Optional

from constants import (
    DEFAULT_DISEASE_MODE,
    DEFAULT_OUTPUT_DIR,
    GRCH38_FASTA,
    resolve_disease_mode,
)
from clinvar_phenotype_matcher import (
    enrich_variants_with_clinvar,
    filter_variants_by_clinvar_disease,
)
from contribution_scorer import ContributionResult, score as score_contribution
from disease_profile import DiseaseProfile
from disease_profile_builder import build_or_load_profile
from disease_reference import get_disease_reference
from disease_space_query import query_disease_space
from gpa_runner import run_gpa_on_filtered_vcf
from hpo_mapper import resolve_disease_query
from liftover import detect_genome_build, liftover_vcf, validate_ref_alleles
from report import generate_report
from tier_filters import denoise_tier3, remove_duplicate_variants
from variant_domain_dive import run_domain_dive_for_variants
from vcf_filter import build_gene_coordinates_cache, check_vcf_completeness, detect_vcf_genotyping

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
    gwas_enabled: bool = True
    literature_enabled: bool = True
    assume_genotyped: bool = False
    assume_not_genotyped: bool = False


def run_disease_risk_pipeline(config: PipelineConfig) -> dict:
    """Execute full Step 1-6 pipeline and return results + report path."""
    output_dir = config.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    disease_mode = resolve_disease_mode(config.disease_query, config.disease_mode)
    logger.info("Disease mode resolved to: %s", disease_mode)

    run_id = f"drq_{config.disease_query.replace(' ', '_').replace('/', '_')[:30]}"
    work_dir = output_dir / run_id
    work_dir.mkdir(parents=True, exist_ok=True)

    # --- Global pipeline checkpoint ---
    result_json = work_dir / "result.json"
    marker = work_dir / ".pipeline_marker"
    if result_json.exists() and marker.exists():
        try:
            input_stat = config.input_vcf.stat()
            stored = json.loads(marker.read_text())
            if (stored.get("input_mtime") == input_stat.st_mtime
                    and stored.get("input_size") == input_stat.st_size
                    and stored.get("disease") == config.disease_query
                    and stored.get("offline") == config.offline
                    and stored.get("spliceai") == config.spliceai):
                logger.info("Pipeline checkpoint hit — returning cached %s", result_json)
                cached = json.loads(result_json.read_text())
                cached["_from_checkpoint"] = True
                return cached
            else:
                logger.info("Pipeline checkpoint stale — re-running")
                marker.unlink(missing_ok=True)
        except Exception as exc:
            logger.warning("Pipeline checkpoint read failed: %s — re-running", exc)
            marker.unlink(missing_ok=True)

    # Step 0: VCF genotyping status detection
    logger.info("Step 0: Detecting VCF genotyping status")
    genotyping_info = detect_vcf_genotyping(config.input_vcf)

    # Apply user overrides
    if config.assume_genotyped:
        genotyping_info["is_genotyped"] = True
        genotyping_info["method"] = "user_override_genotyped"
        genotyping_info["note"] = (
            "User explicitly specified --assume-genotyped. "
            "Missing positions will be treated as REF/REF."
        )
        logger.info("Genotyping status overridden by user: genotyped")
    elif config.assume_not_genotyped:
        genotyping_info["is_genotyped"] = False
        genotyping_info["method"] = "user_override_not_genotyped"
        genotyping_info["note"] = (
            "User explicitly specified --assume-not-genotyped. "
            "Missing positions CANNOT be safely treated as REF/REF."
        )
        logger.info("Genotyping status overridden by user: NOT genotyped")

    is_genotyped = genotyping_info["is_genotyped"]
    logger.info("VCF genotyping status: is_genotyped=%s, method=%s", is_genotyped, genotyping_info["method"])

    if not is_genotyped:
        logger.warning(
            "⚠️ VCF appears NOT genotyped: %s. "
            "Score interpretation may be unreliable — "
            "missing positions may be data gaps, not REF/REF.",
            genotyping_info.get("note", ""),
        )

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

    # Step 2/3: Disease -> HPO and DiseaseProfile
    logger.info("Step 2/3: Mapping disease query to HPO and building DiseaseProfile")
    hpo_result = resolve_disease_query(config.disease_query, config.hpo_id)
    hpo_id = hpo_result.get("hpo_id")
    hpo_name = hpo_result.get("hpo_name")
    if hpo_id is None:
        logger.warning("No HPO mapping found for '%s'; proceeding with OMIM/template only", config.disease_query)

    profile = build_or_load_profile(
        config.disease_query,
        hpo_id=hpo_id,
        refresh=config.refresh_cache,
    )
    profile.to_json(work_dir / "profile.json")
    logger.info(
        "DiseaseProfile loaded: %d genes, %d known variants, %d regulatory regions, mode=%s",
        len(profile.gene_set),
        len(profile.all_known_variants),
        len(profile.regulatory_regions),
        profile.mode,
    )

    # Backward-compatible disease reference cache (for report/legacy fields)
    disease_ref = get_disease_reference(config.disease_query, refresh=config.refresh_cache)
    # Sync literature count with the dynamic profile literature
    if isinstance(disease_ref, dict) and "metadata" in disease_ref and "counts" in disease_ref["metadata"]:
        disease_ref["metadata"]["counts"]["literature_entries"] = len(profile.key_literature or [])

    if not profile.gene_set:
        # No genes found; produce empty report
        report_path = work_dir / f"{config.disease_query.replace(' ', '_')}-{datetime.now().strftime('%Y%m%d')}-report.md"
        empty_gpa = {
            "tier1_variants": [],
            "tier2_variants": [],
            "tier3_variants": [],
            "multi_hit": [],
            "summary": {},
        }
        empty_contribution = ContributionResult(
            overall_level="uncertain",
            overall_score=None,
        )
        gene_set_result = {
            "merged_genes": [],
            "sources": {},
            "total": 0,
            "disease_reference": disease_ref,
        }
        generate_report(
            disease_name=config.disease_query,
            hpo_id=hpo_id,
            hpo_name=hpo_name,
            sex=config.sex,
            age=config.age,
            gene_set_result=gene_set_result,
            gpa_result=empty_gpa,
            score_result=empty_contribution.to_dict(),
            output_path=report_path,
            apoe_result=None,
            vcf_qc=vcf_qc,
            disease_mode=disease_mode,
            is_genotyped=is_genotyped,
            genotyping_info=genotyping_info,
        )
        return {
            "success": True,
            "hpo_id": hpo_id,
            "hpo_name": hpo_name,
            "disease_profile": profile.to_dict(),
            "gene_set": gene_set_result,
            "gpa": empty_gpa,
            "contribution": empty_contribution.to_dict(),
            "score": empty_contribution.to_dict(),
            "apoe": None,
            "report_path": report_path,
            "work_dir": work_dir,
            "validation": validation,
            "note": "No disease-associated genes found in local databases.",
        }

    # Step 4: Unified disease-space query
    logger.info("Step 4: Querying unified disease space")
    gene_coords = build_gene_coordinates_cache()
    space_result = query_disease_space(
        normalized_vcf,
        profile,
        gene_coords,
        work_dir,
    )
    filtered_vcf = space_result["disease_space_vcf"]
    known_genotypes = space_result["known_variant_genotypes"]
    logger.info(
        "Disease space variants: %d, known variants queried: %d",
        space_result["variant_count"], len(known_genotypes),
    )

    # Build legacy gene_set_result for report compatibility
    hpo_genes = set(hpo_result.get("genes", [])) if hpo_id else set()
    builtin_genes = {g.gene for g in profile.gene_set}
    omim_only = set()
    # OMIM-derived genes are those not from built-in template and not from HPO
    # (heuristic; exact source depends on builder path)
    hpo_only = hpo_genes - builtin_genes
    extra_genes = set()
    for g in profile.gene_set:
        if g.evidence == "hpo_omim":
            extra_genes.add(g.gene)
        if g.evidence == "omim":
            omim_only.add(g.gene)
    omim_only = omim_only - builtin_genes - hpo_genes
    extra_genes = extra_genes - builtin_genes - hpo_genes

    gene_set_result = {
        "merged_genes": sorted(profile.all_genes),
        "sources": {
            "omim": len(omim_only),
            "hpo": len(hpo_only),
            "extra": len(extra_genes),
            "builtin": len(builtin_genes),
        },
        "total": len(profile.all_genes),
        "disease_reference": disease_ref,
    }

    # Step 5A/B: GPA tier classification on disease-space VCF
    logger.info("Step 5: Running GPA tier classification")
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
            "hpo_name": hpo_name,
            "disease_profile": profile.to_dict(),
            "gene_set": gene_set_result,
            "work_dir": work_dir,
        }

    # Step 5C: Post-GPA filtering
    logger.info("Step 5C: Post-GPA ClinVar phenotype matching and Tier 3 denoising")
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
        core_genes=profile.core_genes,
        gwas_loci_genes=profile.gwas_loci_genes,
    )

    # Step 5C-ter: Enrich all tiered variants with ClinVar annotations.
    # This is informational only: P/LP is highlighted, but VUS/conflicting/etc.
    # are not used to filter or downweight variants.
    logger.info("Step 5C-ter: Enriching tiered variants with ClinVar annotations")
    all_tiered = (
        gpa_result.get("tier1_variants", [])
        + gpa_result.get("tier2_variants", [])
        + gpa_result.get("tier3_variants", [])
    )
    enriched = enrich_variants_with_clinvar(all_tiered, config.disease_query)
    # Split back into tiers preserving order
    t1_n = len(gpa_result.get("tier1_variants", []))
    t2_n = len(gpa_result.get("tier2_variants", []))
    gpa_result["tier1_variants"] = enriched[:t1_n]
    gpa_result["tier2_variants"] = enriched[t1_n:t1_n + t2_n]
    gpa_result["tier3_variants"] = enriched[t1_n + t2_n:]
    logger.info("ClinVar enrichment completed for %d variants", len(enriched))

    # Step 5C-bis: Targeted domain-dive for Tier 2/3 variants in core genes
    logger.info("Step 5C-bis: Running targeted domain-dive on Tier 2/3 candidates")
    tier2_and_3 = (
        gpa_result.get("tier2_variants", [])
        + gpa_result.get("tier3_variants", [])
    )
    # Pass key_regions to domain dive via existing mechanism if available
    domain_dive_candidates = run_domain_dive_for_variants(
        tier2_and_3,
        disease_name=config.disease_query,
        min_tier=2,
        max_tier=3,
    )
    gpa_result["domain_dive_candidates"] = domain_dive_candidates
    logger.info("Domain-dive candidates: %d", len(domain_dive_candidates))

    # Rebuild multi-hit from filtered tiers
    from collections import Counter
    gene_counts = Counter(v.get("gene") or v.get("GENE") for v in (
        gpa_result.get("tier1_variants", [])
        + gpa_result.get("tier2_variants", [])
        + gpa_result.get("tier3_variants", [])
    ) if (v.get("gene") or v.get("GENE")))
    gpa_result["multi_hit"] = [
        {"gene": gene, "variant_count": count}
        for gene, count in gene_counts.items()
        if count > 1
    ]

    # Step 6: Contribution scoring
    logger.info("Step 6: Calculating disease contribution score")
    contribution = score_contribution(
        profile=profile,
        tiered_variants={
            "tier1_variants": gpa_result.get("tier1_variants", []),
            "tier2_variants": gpa_result.get("tier2_variants", []),
            "tier3_variants": gpa_result.get("tier3_variants", []),
        },
        known_genotypes=known_genotypes,
        vcf_qc=vcf_qc,
    )
    score_result = contribution.to_dict()
    # Backward-compatible report keys
    score_result["total_score"] = round((contribution.overall_score or 0.0) * 100)
    score_result["risk_level"] = contribution.overall_level
    score_result["contribution_level"] = contribution.overall_level
    score_result["risk_meaning"] = contribution.overall_level
    score_result["contribution_meaning"] = contribution.overall_level
    score_result["components"] = {
        "tier1_score": round(sum(x["contribution"] for x in contribution.mendelian_high) * 100),
        "tier2_score": round(sum(x["contribution"] for x in contribution.mendelian_mod) * 100),
        "monogenic_score": round(sum(x["contribution"] for x in contribution.mendelian_high) * 100),
        "rare_functional_score": round(sum(x["contribution"] for x in contribution.mendelian_mod) * 100),
        "gwas_score": round(abs(contribution.gwas_prs.get("score", 0.0)) * 100),
        "gwas_hits": contribution.gwas_prs.get("variants", []),
        "literature_score": 0,
        "rarity_score": 0.0,
    }
    logger.info(
        "Contribution score: %.3f, level: %s",
        score_result.get("overall_score") or 0,
        score_result["overall_level"],
    )

    # Apply literature/GWAS flags to variants for report
    for v in gpa_result.get("tier1_variants", []) + gpa_result.get("tier2_variants", []):
        flags = []
        gene = v.get("gene") or v.get("GENE") or ""
        if gene in profile.gwas_loci_genes:
            flags.append("GWAS_LOCUS")
        v["_gene_contribution"] = profile.gene_contribution_map.get(gene)
        v["_gene_penetrance"] = profile.gene_penetrance_map.get(gene)
        v["_drq_flags"] = flags

    # Build literature summary from profile key_literature
    literature_entries = profile.key_literature or []
    literature_genes = set()
    for entry in literature_entries:
        for gene in entry.get("genes", []):
            literature_genes.add(gene.upper())
    literature_summary = {
        "variant_hits": [],
        "gene_hits": sorted(literature_genes),
        "core_genes_covered": sorted(literature_genes & profile.core_genes),
        "total_entries": len(literature_entries),
        "entries": literature_entries,
    }

    # Generate report
    report_path = work_dir / f"{config.disease_query.replace(' ', '_')}-{datetime.now().strftime('%Y%m%d')}-report.md"
    legacy_gwas_lead_snps = _known_genotypes_to_legacy_gwas(known_genotypes)

    # Build disease space stats for report
    known_found = sum(1 for kg in known_genotypes if kg.dosage > 0 or not kg.inferred_ref_ref)
    known_inferred = sum(1 for kg in known_genotypes if kg.inferred_ref_ref)
    disease_space = {
        "total_variants": vcf_qc.get("total_variants", 0) if vcf_qc else 0,
        "analyzed_variants": space_result["variant_count"],
        "known_variants_queried": len(known_genotypes),
        "known_found": known_found,
        "known_inferred": known_inferred,
    }

    generate_report(
        disease_name=config.disease_query,
        hpo_id=hpo_id,
        hpo_name=hpo_name,
        sex=config.sex,
        age=config.age,
        gene_set_result=gene_set_result,
        gpa_result=gpa_result,
        score_result=score_result,
        apoe_result=None,
        gwas_summary={"hit_count": 0, "hit_genes": [], "top_variants": []},
        literature_summary=literature_summary,
        disease_reference=disease_ref,
        gwas_lead_snps=legacy_gwas_lead_snps,
        vcf_qc=vcf_qc,
        output_path=report_path,
        disease_mode=disease_mode,
        domain_dive_candidates=domain_dive_candidates,
        disease_space=disease_space,
        is_genotyped=is_genotyped,
        genotyping_info=genotyping_info,
    )

    # Save structured JSON
    json_path = work_dir / "result.json"
    result = {
        "success": True,
        "hpo_id": hpo_id,
        "hpo_name": hpo_name,
        "is_genotyped": is_genotyped,
        "genotyping_info": genotyping_info,
        "input_vcf": str(config.input_vcf),
        "normalized_vcf": str(normalized_vcf),
        "filtered_vcf": str(filtered_vcf),
        "genome_build": build,
        "validation": validation,
        "apoe": None,
        "disease_mode": disease_mode,
        "disease_profile": profile.to_dict(),
        "gene_set": gene_set_result,
        "filter_stats": {
            "input_variants": space_result["variant_count"],
            "output_variants": space_result["variant_count"],
            "gene_count": len(profile.all_genes),
            "known_variant_count": len(known_genotypes),
        },
        "gpa": gpa_result,
        "contribution": contribution.to_dict(),
        "score": score_result,
        "domain_dive_candidates": domain_dive_candidates,
        "known_variant_genotypes": [g.to_dict() for g in known_genotypes],
        "gwas_lead_snps": legacy_gwas_lead_snps,
        "gwas_summary": {"hit_count": 0, "hit_genes": [], "top_variants": []},
        "literature_summary": literature_summary,
        "disease_reference": disease_ref,
        "vcf_qc": vcf_qc,
        "report_path": str(report_path),
        "work_dir": str(work_dir),
    }
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, ensure_ascii=False, default=str)

    # Write pipeline checkpoint marker
    try:
        input_stat = config.input_vcf.stat()
        marker.write_text(
            json.dumps(
                {
                    "input_mtime": input_stat.st_mtime,
                    "input_size": input_stat.st_size,
                    "disease": config.disease_query,
                    "offline": config.offline,
                    "spliceai": config.spliceai,
                }
            )
        )
    except Exception as exc:
        logger.warning("Pipeline marker write failed (non-fatal): %s", exc)

    return result


def _known_genotypes_to_legacy_gwas(known_genotypes: list) -> list[dict]:
    """Convert KnownVariantGenotype objects to legacy gwas_lead_snps format."""
    results = []
    for kg in known_genotypes:
        v = kg.variant
        results.append({
            "rsid": v.rsid,
            "gene": v.gene,
            "chrom": kg.chrom,
            "pos": kg.pos,
            "ref": kg.ref,
            "alt": kg.alt,
            "effect_allele": v.effect_allele,
            "beta": v.beta,
            "or": v.or_value,
            "contribution_score": v.contribution_score,
            "sample_gt": {
                "chrom": kg.chrom,
                "pos": kg.pos,
                "ref": kg.ref,
                "alt": kg.alt,
                "gt": kg.gt,
                "filter": kg.filter_status,
                "format": kg.sample_format,
                "sample": kg.sample_values,
            },
            "note": v.note,
        })
    return results
