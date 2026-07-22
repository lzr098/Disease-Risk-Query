#!/usr/bin/env python
"""Quick template-only test for psychiatric disorders composite template.

Only queries the 196 known variants (GWAS lead SNPs + PRS SNPs) from the VCF,
skipping the full disease-space extraction (71K variants) and GPA tier
classification.  This validates:
  - Template loading (184 genes, 196 known variants)
  - Known variant genotype query (bcftools)
  - Contribution scoring with sub_disease breakdown
  - Report generation (including psychiatric disclaimer + sub-disease table)
"""

import json
import logging
import sys
from datetime import datetime
from pathlib import Path

# Ensure scripts/ is on the path
SCRIPTS_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SCRIPTS_DIR))

from disease_profile_builder import build_or_load_profile
from disease_space_query import query_known_variants
from contribution_scorer import ContributionResult, score as score_contribution
from report import generate_report
from vcf_filter import build_gene_coordinates_cache
from hpo_mapper import resolve_disease_query

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger(__name__)

VCF_PATH = Path("/Users/zhaorongli/Documents/Documents - zhaorong\u2019s MacBook Air/\u6d3b\u52a8\u8d44\u6599/F57J018690-WGS.final.vcf.gz")
OUTPUT_DIR = Path("/tmp/drq_psych_test_known_only")
DISEASE = "psychiatric disorders"
SEX = "female"
AGE = 10


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Step 1: Resolve HPO
    logger.info("Step 1: Resolving HPO for '%s'", DISEASE)
    hpo_result = resolve_disease_query(DISEASE)
    hpo_id = hpo_result.get("hpo_id", "")
    hpo_name = hpo_result.get("hpo_name", "")
    logger.info("HPO: %s (%s)", hpo_id, hpo_name)

    # Step 2: Build DiseaseProfile
    logger.info("Step 2: Building DiseaseProfile")
    profile = build_or_load_profile(DISEASE, hpo_id=hpo_id)
    logger.info("Profile: %d genes, %d known variants, mode=%s",
                len(profile.gene_set), len(profile.all_known_variants), profile.mode)

    # Step 3: Query ONLY known variants (skip disease-space VCF extraction)
    logger.info("Step 3: Querying %d known variant positions", len(profile.all_known_variants))
    known_genotypes = query_known_variants(VCF_PATH, profile.all_known_variants, sex=SEX)
    found = sum(1 for kg in known_genotypes if kg.dosage > 0)
    inferred = sum(1 for kg in known_genotypes if kg.inferred_ref_ref)
    logger.info("Known variants: %d found, %d inferred ref/ref, %d total",
                found, inferred, len(known_genotypes))

    # Print sub-disease breakdown
    from collections import Counter
    sub_found = Counter()
    sub_inferred = Counter()
    for kg in known_genotypes:
        sd = kg.variant.sub_disease or "unknown"
        if kg.dosage > 0:
            sub_found[sd] += 1
        if kg.inferred_ref_ref:
            sub_inferred[sd] += 1
    logger.info("Sub-disease breakdown (found / inferred):")
    for sd in sorted(sub_found | sub_inferred):
        logger.info("  %s: %d found / %d inferred", sd, sub_found.get(sd, 0), sub_inferred.get(sd, 0))

    # Step 4: Contribution scoring (empty GPA tiers - template-only test)
    logger.info("Step 4: Running contribution scorer")
    vcf_qc = {
        "checked": True,
        "anchor_snps_checked": len(profile.all_known_variants),
        "anchor_snps_present": found,
        "presence_rate": found / max(len(profile.all_known_variants), 1),
        "common_variants_filtered": False,
        "total_variants": 5180075,
        "threshold": 0.5,
        "note": "pass",
    }
    contribution = score_contribution(
        profile=profile,
        tiered_variants={
            "tier1_variants": [],
            "tier2_variants": [],
            "tier3_variants": [],
        },
        known_genotypes=known_genotypes,
        vcf_qc=vcf_qc,
    )
    score_result = contribution.to_dict()
    score_result["total_score"] = round((contribution.overall_score or 0.0) * 100)
    score_result["risk_level"] = contribution.overall_level
    score_result["contribution_level"] = contribution.overall_level
    score_result["risk_meaning"] = contribution.overall_level
    score_result["contribution_meaning"] = contribution.overall_level
    score_result["components"] = {
        "tier1_score": 0,
        "tier2_score": 0,
        "monogenic_score": 0,
        "rare_functional_score": 0,
        "gwas_score": round(abs(contribution.gwas_prs.get("score", 0.0)) * 100),
        "gwas_hits": contribution.gwas_prs.get("variants", []),
        "literature_score": 0,
        "rarity_score": 0.0,
    }
    logger.info("Overall score: %.3f, level: %s",
                score_result.get("overall_score") or 0, score_result["overall_level"])

    # Print sub-disease scores
    if hasattr(contribution, 'sub_disease_scores') and contribution.sub_disease_scores:
        logger.info("Sub-disease scores:")
        for sd, sc in sorted(contribution.sub_disease_scores.items()):
            logger.info("  %s: score=%.4f, level=%s", sd, sc.get("score", 0), sc.get("level", "N/A"))

    # Step 5: Build gene_set_result
    builtin_genes = {g.gene for g in profile.gene_set}
    gene_set_result = {
        "merged_genes": sorted(profile.all_genes),
        "sources": {"omim": 0, "hpo": 0, "extra": 0, "builtin": len(builtin_genes)},
        "total": len(profile.all_genes),
        "disease_reference": None,
    }

    # Step 6: Build legacy known genotypes for report
    from pipeline import _known_genotypes_to_legacy_gwas
    legacy_gwas = _known_genotypes_to_legacy_gwas(known_genotypes)

    # Build disease space stats
    disease_space = {
        "total_variants": 5180075,
        "analyzed_variants": 0,  # Skipped disease-space extraction
        "known_variants_queried": len(known_genotypes),
        "known_found": found,
        "known_inferred": inferred,
    }

    # Build literature summary
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

    # Step 7: Generate report
    report_path = OUTPUT_DIR / f"psychiatric_disorders-{datetime.now().strftime('%Y%m%d_%H%M%S')}-report.md"
    logger.info("Step 5: Generating report -> %s", report_path)

    gpa_result = {
        "success": True,
        "tier1_variants": [],
        "tier2_variants": [],
        "tier3_variants": [],
        "multi_hit": [],
        "domain_dive_candidates": [],
    }

    generate_report(
        disease_name=DISEASE,
        hpo_id=hpo_id,
        hpo_name=hpo_name,
        sex=SEX,
        age=AGE,
        gene_set_result=gene_set_result,
        gpa_result=gpa_result,
        score_result=score_result,
        apoe_result=None,
        gwas_summary={"hit_count": 0, "hit_genes": [], "top_variants": []},
        literature_summary=literature_summary,
        disease_reference=None,
        gwas_lead_snps=legacy_gwas,
        vcf_qc=vcf_qc,
        output_path=report_path,
        disease_mode="complex",
        domain_dive_candidates=[],
        disease_space=disease_space,
        is_genotyped=True,
        genotyping_info={"method": "auto", "note": "Template-only test, disease-space skipped"},
    )

    # Save structured JSON
    json_path = OUTPUT_DIR / "result.json"
    result = {
        "success": True,
        "disease": DISEASE,
        "profile_genes": len(profile.gene_set),
        "known_variants": len(profile.all_known_variants),
        "known_found": found,
        "known_inferred": inferred,
        "sub_disease_found": dict(sub_found),
        "sub_disease_inferred": dict(sub_inferred),
        "overall_score": score_result.get("overall_score"),
        "overall_level": score_result["overall_level"],
        "sub_disease_scores": getattr(contribution, 'sub_disease_scores', {}),
        "report_path": str(report_path),
    }
    json_path.write_text(json.dumps(result, indent=2, ensure_ascii=False))
    logger.info("Done! Report: %s", report_path)
    logger.info("JSON: %s", json_path)
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
