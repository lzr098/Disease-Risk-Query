"""Validate built-in disease templates for structural and data-quality issues.

Enhanced v2: Added EAF coverage, KPV gene coverage, gene annotation completeness,
ClinGen validity coverage, key_regions coverage, CS cross-disease consistency,
literature freshness, and irrelevant gene detection.
"""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from constants import DISEASE_BUILTIN_REFS
from disease_profile import DiseaseProfile, VariantWeight


def _variant_key(v: dict) -> str:
    return f"{v.get('chrom')}:{v.get('pos')}:{v.get('ref', '')}>{v.get('alt', '')}"


def validate() -> list[str]:
    issues: list[str] = []
    warnings: list[str] = []

    for disease, ref in DISEASE_BUILTIN_REFS.items():
        # ---- 1) Required fields ----
        required = [
            "aliases", "mode", "gene_set", "known_pathogenic_variants",
            "gwas_lead_snps", "prs_variants", "prs_variants_high",
            "regulatory_regions", "key_regions", "key_literature",
        ]
        for field in required:
            if field not in ref:
                issues.append(f"{disease}: missing required field '{field}'")

        # ---- 2) Mode validity ----
        mode = ref.get("mode")
        if mode not in ("mendelian", "complex", "mixed"):
            issues.append(f"{disease}: invalid mode '{mode}'")

        # ---- 3) Aliases ----
        aliases = ref.get("aliases", [])
        if len(aliases) < 3:
            warnings.append(f"{disease}: only {len(aliases)} aliases (recommended ≥3)")

        # ---- 4) Duplicate genes ----
        gene_set = ref.get("gene_set", [])
        genes = [g.get("gene") for g in gene_set]
        dup_genes = {g for g, c in Counter(genes).items() if c > 1}
        if dup_genes:
            issues.append(f"{disease}: duplicate genes: {sorted(dup_genes)}")

        # ---- 5) Gene field completeness ----
        for g in gene_set:
            if not g.get("gene"):
                issues.append(f"{disease}: gene entry missing gene symbol")
            if not g.get("tier"):
                issues.append(f"{disease}: gene '{g.get('gene')}' missing tier")
            if g.get("tier") not in ("mendelian_high", "mendelian_mod",
                                      "strong_gwas", "gwas", "dosage"):
                issues.append(f"{disease}: gene '{g.get('gene')}' invalid tier '{g.get('tier')}'")
            cs = g.get("contribution_score", 0)
            if cs and (cs < 0.05 or cs > 1.1):
                issues.append(f"{disease}: gene '{g.get('gene')}' CS={cs} out of range [0.05, 1.0]")

        # ---- 6) Gene annotation coverage ----
        mendelian_genes = [g for g in gene_set if g.get("tier") in ("mendelian_high", "mendelian_mod")]
        for g in mendelian_genes:
            if not g.get("phenotype_assoc"):
                warnings.append(f"{disease}: '{g.get('gene')}' (mendelian) missing phenotype_assoc")
            if not g.get("key_domains"):
                warnings.append(f"{disease}: '{g.get('gene')}' (mendelian) missing key_domains")

        clingen_annotated = sum(1 for g in mendelian_genes if g.get("clingen_validity"))
        if mendelian_genes and clingen_annotated / len(mendelian_genes) < 0.5:
            warnings.append(
                f"{disease}: only {clingen_annotated}/{len(mendelian_genes)} "
                f"mendelian genes have clingen_validity"
            )

        # ---- 7) Irrelevant gene detection ----
        weak_keywords = ["继发性", "间接", "可能相关", "secondary", "indirect"]
        for g in gene_set:
            note = g.get("note", "")
            for kw in weak_keywords:
                if kw in note:
                    warnings.append(
                        f"{disease}: '{g.get('gene')}' note contains '{kw}' — "
                        f"consider removing or downgrading tier"
                    )
                    break

        # ---- 8) Variant duplicate detection ----
        all_variants = (
            list(ref.get("known_pathogenic_variants", []))
            + list(ref.get("gwas_lead_snps", []))
            + list(ref.get("prs_variants", []))
            + list(ref.get("prs_variants_high", []))
        )
        variant_ids = [
            (v.get("rsid"), v.get("chrom"), v.get("pos"), v.get("ref"), v.get("alt"))
            for v in all_variants
        ]
        dup_variant_ids = {vid for vid, c in Counter(variant_ids).items() if c > 1}
        if dup_variant_ids:
            issues.append(f"{disease}: duplicate variant entries: {sorted(dup_variant_ids)}")

        # ---- 9) Variant field completeness ----
        for list_name in ("known_pathogenic_variants", "gwas_lead_snps",
                           "prs_variants", "prs_variants_high"):
            for v in ref.get(list_name, []):
                for field in ("chrom", "pos"):
                    if field not in v or v.get(field) is None:
                        issues.append(
                            f"{disease}.{list_name}: variant missing {field}: {v.get('rsid')}"
                        )
                if "or" in v and "or_value" not in v:
                    issues.append(
                        f"{disease}.{list_name}: legacy 'or' key: {v.get('rsid')}"
                    )
                try:
                    VariantWeight.from_dict(v)
                except Exception as exc:
                    issues.append(
                        f"{disease}.{list_name}: VariantWeight failed for {v.get('rsid')}: {exc}"
                    )

        # ---- 10) GWAS SNP EAF coverage ----
        gwas_snps = ref.get("gwas_lead_snps", [])
        if gwas_snps:
            eur_filled = sum(1 for s in gwas_snps if s.get("eaf_eur") is not None)
            eas_filled = sum(1 for s in gwas_snps if s.get("eaf_eas") is not None)
            eur_pct = eur_filled / len(gwas_snps) * 100
            eas_pct = eas_filled / len(gwas_snps) * 100
            if eur_pct < 70:
                issues.append(
                    f"{disease}: EUR EAF coverage {eur_pct:.0f}% < 70% "
                    f"({eur_filled}/{len(gwas_snps)})"
                )
            if eas_pct < 60:
                issues.append(
                    f"{disease}: EAS EAF coverage {eas_pct:.0f}% < 70% "
                    f"({eas_filled}/{len(gwas_snps)})"
                )

            # Effect size completeness
            no_effect = sum(1 for s in gwas_snps
                           if s.get("or_value") is None and s.get("beta") is None)
            if no_effect > 0:
                issues.append(
                    f"{disease}: {no_effect}/{len(gwas_snps)} GWAS SNPs "
                    f"missing both or_value and beta"
                )

        # ---- 11) GWAS SNP CS cross-reference ----
        cs_vals = [s.get("contribution_score", 0) for s in gwas_snps if s.get("contribution_score")]
        if cs_vals:
            mean_cs = sum(cs_vals) / len(cs_vals)
            if mean_cs < 0.10 or mean_cs > 0.35:
                issues.append(
                    f"{disease}: GWAS CS mean={mean_cs:.3f} outside [0.10, 0.35]"
                )

        # ---- 12) KPV gene coverage ----
        kpv_genes = set(v.get("gene", "") for v in ref.get("known_pathogenic_variants", []))
        mh_genes = [g["gene"] for g in gene_set if g.get("tier") == "mendelian_high"]
        mm_genes = [g["gene"] for g in gene_set if g.get("tier") == "mendelian_mod"]
        for gene in mh_genes:
            if gene not in kpv_genes:
                warnings.append(f"{disease}: mendelian_high gene '{gene}' has no KPV")
        mm_no_kpv = sum(1 for g in mm_genes if g not in kpv_genes)
        if mm_no_kpv > len(mm_genes) * 0.5 and mm_genes:
            warnings.append(
                f"{disease}: {mm_no_kpv}/{len(mm_genes)} mendelian_mod genes without KPV"
            )

        # ---- 13) key_regions coverage ----
        kr = ref.get("key_regions", {})
        if isinstance(kr, dict):
            if len(kr) < 2:
                warnings.append(f"{disease}: key_regions covers only {len(kr)} genes (recommended ≥2)")
        else:
            warnings.append(f"{disease}: key_regions is not a dict")

        # ---- 14) PRS coverage ----
        prs_high = ref.get("prs_variants_high", [])
        if mode in ("complex", "mixed") and len(prs_high) < 3:
            warnings.append(
                f"{disease}: prs_variants_high has {len(prs_high)} entries "
                f"(recommended ≥3 for {mode} mode)"
            )

        # ---- 15) Literature freshness ----
        lit = ref.get("key_literature", [])
        if len(lit) < 3:
            warnings.append(f"{disease}: only {len(lit)} literature entries (recommended ≥3)")
        recent = sum(1 for l in lit if int(str(l.get("pmid", "0")[:4]) or "0") >= 2020
                    or "2020" in l.get("note", "") or "2021" in l.get("note", "")
                    or "2022" in l.get("note", ""))
        if recent == 0 and lit:
            warnings.append(f"{disease}: no literature from 2020+")

        # ---- 16) regulatory_regions ----
        for r in ref.get("regulatory_regions", []):
            for field in ("chrom", "start", "end"):
                if field not in r:
                    issues.append(f"{disease}.regulatory_regions: missing {field}")

        # ---- 17) Full profile conversion ----
        try:
            profile = DiseaseProfile.from_dict({"canonical_name": disease, **ref})
            _ = profile.to_dict()
        except Exception as exc:
            issues.append(f"{disease}: DiseaseProfile conversion failed: {exc}")

    # ---- 18) Cross-disease CS consistency ----
    tier_cs: dict[str, list[float]] = {}
    for disease, ref in DISEASE_BUILTIN_REFS.items():
        for g in ref.get("gene_set", []):
            tier = g.get("tier", "")
            cs = g.get("contribution_score", 0)
            if tier and cs:
                tier_cs.setdefault(tier, []).append(cs)

    for tier, cs_list in tier_cs.items():
        if len(cs_list) < 3:
            continue
        mean_cs = sum(cs_list) / len(cs_list)
        variance = sum((x - mean_cs) ** 2 for x in cs_list) / len(cs_list)
        std = variance ** 0.5
        if std > 0.15:
            warnings.append(
                f"cross-disease: tier '{tier}' CS std={std:.3f} > 0.15 "
                f"(mean={mean_cs:.3f}, n={len(cs_list)})"
            )

    # Report: issues first, then warnings
    result = []
    if issues:
        result.append(f"=== {len(issues)} ERROR(S) ===")
        for i in issues:
            result.append(f"  [ERR] {i}")
    if warnings:
        result.append(f"=== {len(warnings)} WARNING(S) ===")
        for w in warnings:
            result.append(f"  [WARN] {w}")
    if not issues and not warnings:
        result.append("All built-in disease templates passed validation.")

    return result


def main() -> int:
    lines = validate()
    for line in lines:
        print(line)
    has_errors = any("[ERR]" in line for line in lines)
    return 1 if has_errors else 0


if __name__ == "__main__":
    sys.exit(main())
