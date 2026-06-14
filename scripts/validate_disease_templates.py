"""Validate built-in disease templates for structural and data-quality issues."""

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

    for disease, ref in DISEASE_BUILTIN_REFS.items():
        # Required fields
        required = [
            "mode", "gene_set", "known_pathogenic_variants", "gwas_lead_snps",
            "prs_variants", "regulatory_regions", "key_regions", "key_literature",
        ]
        for field in required:
            if field not in ref:
                issues.append(f"{disease}: missing required field '{field}'")

        # Mode validity
        mode = ref.get("mode")
        if mode not in ("mendelian", "complex", "mixed"):
            issues.append(f"{disease}: invalid mode '{mode}'")

        # Duplicate genes
        gene_set = ref.get("gene_set", [])
        genes = [g.get("gene") for g in gene_set]
        dup_genes = {g for g, c in Counter(genes).items() if c > 1}
        if dup_genes:
            issues.append(f"{disease}: duplicate genes: {sorted(dup_genes)}")

        # Genes with missing/empty tier
        for g in gene_set:
            if not g.get("gene"):
                issues.append(f"{disease}: gene entry missing gene symbol")
            if not g.get("tier"):
                issues.append(f"{disease}: gene '{g.get('gene')}' missing tier")

        # Duplicate variant entries across all variant lists
        all_variants = (
            list(ref.get("known_pathogenic_variants", []))
            + list(ref.get("gwas_lead_snps", []))
            + list(ref.get("prs_variants", []))
        )
        variant_ids = [
            (v.get("rsid"), v.get("chrom"), v.get("pos"), v.get("ref"), v.get("alt"))
            for v in all_variants
        ]
        dup_variant_ids = {vid for vid, c in Counter(variant_ids).items() if c > 1}
        if dup_variant_ids:
            issues.append(f"{disease}: duplicate variant entries: {sorted(dup_variant_ids)}")

        # Variant field completeness and VariantWeight compatibility
        for list_name in ("known_pathogenic_variants", "gwas_lead_snps", "prs_variants"):
            for v in ref.get(list_name, []):
                for field in ("chrom", "pos"):
                    if field not in v or v.get(field) is None:
                        issues.append(f"{disease}.{list_name}: variant missing {field}: {v.get('rsid')}")
                # Check for legacy 'or' key instead of 'or_value'
                if "or" in v and "or_value" not in v:
                    issues.append(f"{disease}.{list_name}: variant uses legacy 'or' key: {v.get('rsid')}")
                # Try converting to VariantWeight
                try:
                    VariantWeight.from_dict(v)
                except Exception as exc:
                    issues.append(f"{disease}.{list_name}: VariantWeight.from_dict failed for {v.get('rsid')}: {exc}")

        # regulatory_regions
        for r in ref.get("regulatory_regions", []):
            for field in ("chrom", "start", "end"):
                if field not in r:
                    issues.append(f"{disease}.regulatory_regions: missing {field}")

        # Try full profile conversion
        try:
            profile = DiseaseProfile.from_dict({
                "canonical_name": disease,
                **ref,
            })
            _ = profile.to_dict()
        except Exception as exc:
            issues.append(f"{disease}: DiseaseProfile conversion failed: {exc}")

    return issues


def main() -> int:
    issues = validate()
    if issues:
        print(f"Found {len(issues)} issue(s):")
        for issue in issues:
            print(f"  - {issue}")
        return 1
    print("All built-in disease templates passed validation.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
