"""Targeted variant domain-dive analysis.

Given a Tier 3 (or borderline) variant, check whether it falls into a
 disease-specific key functional region or near a critical residue.
This provides structured evidence for manual Tier 2 upgrade decisions.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Optional

from constants import DISEASE_BUILTIN_REFS, resolve_builtin_disease_key


def parse_protein_position(hgvsp: Optional[str]) -> Optional[int]:
    """Extract the protein position from an HGVSp string.

    Examples:
        ENSP00000368727.3:p.Asn1087Ser -> 1087
        p.Gly920Trp -> 920
        p.Arg881Met -> 881
    """
    if not hgvsp:
        return None
    # Keep only the part after the last ':' if present
    if ":" in hgvsp:
        hgvsp = hgvsp.split(":")[-1]
    # Match position (digits) between two letters, e.g. Asn1087Ser
    match = re.search(r"[A-Za-z]+(\d+)[A-Za-z*=?]+", hgvsp)
    if match:
        return int(match.group(1))
    return None


def _parse_residue_range(range_str: str) -> tuple[int, int]:
    """Parse 'start-end' residue range."""
    start, end = range_str.split("-")
    return int(start), int(end)


def _parse_residue(residue_str: str) -> tuple[str, int]:
    """Parse 'Glu803' -> ('Glu', 803)."""
    match = re.match(r"([A-Za-z]+)(\d+)", residue_str)
    if not match:
        raise ValueError(f"Cannot parse residue: {residue_str}")
    return match.group(1), int(match.group(2))


def get_disease_key_regions(disease_name: str) -> dict[str, dict]:
    """Return key_regions for a disease from built-in templates."""
    key = resolve_builtin_disease_key(disease_name)
    if not key:
        return {}
    return DISEASE_BUILTIN_REFS[key].get("key_regions", {})


def assess_variant_in_key_regions(
    gene: str,
    protein_position: int,
    disease_name: str,
    hgvsp: Optional[str] = None,
) -> dict[str, Any]:
    """Assess whether a variant hits disease-specific key regions.

    Returns a dict with:
      - gene_has_key_regions: bool
      - in_key_region: bool
      - matched_regions: list of matched region dicts
      - nearest_critical_residues: list of nearby critical residues with distances
      - upgrade_recommendation: "tier2_candidate" | "monitor" | "no_evidence"
      - reasoning: str
    """
    key_regions = get_disease_key_regions(disease_name)
    gene_regions = key_regions.get(gene)

    result: dict[str, Any] = {
        "gene": gene,
        "protein_position": protein_position,
        "hgvsp": hgvsp,
        "disease": disease_name,
        "gene_has_key_regions": gene_regions is not None,
        "in_key_region": False,
        "matched_regions": [],
        "nearest_critical_residues": [],
        "upgrade_recommendation": "no_evidence",
        "reasoning": "",
    }

    if not gene_regions:
        result["reasoning"] = f"疾病「{disease_name}」中未定义 {gene} 的关键区域。"
        return result

    # Check region matches
    for region in gene_regions.get("regions", []):
        start, end = _parse_residue_range(region["residues"])
        if start <= protein_position <= end:
            result["in_key_region"] = True
            result["matched_regions"].append(region)

    # Find nearest critical residues
    nearby: list[dict] = []
    for cr in gene_regions.get("critical_residues", []):
        _, pos = _parse_residue(cr["residue"])
        distance = abs(protein_position - pos)
        if distance <= 20:  # within 20 aa window
            nearby.append({
                "residue": cr["residue"],
                "distance": distance,
                "note": cr.get("note", ""),
            })
    nearby.sort(key=lambda x: x["distance"])
    result["nearest_critical_residues"] = nearby

    # Upgrade recommendation logic
    if result["in_key_region"] and nearby:
        result["upgrade_recommendation"] = "tier2_candidate"
        result["reasoning"] = (
            f"{gene} p.{protein_position} 位于疾病关键区域 "
            f"（{', '.join(r['name'] for r in result['matched_regions'])}），"
            f"且邻近关键残基 "
            f"（{', '.join(f'{n['residue']}({n['distance']}aa)' for n in nearby[:3])}），"
            f"建议人工复核是否升级为 Tier 2。"
        )
    elif result["in_key_region"]:
        result["upgrade_recommendation"] = "monitor"
        result["reasoning"] = (
            f"{gene} p.{protein_position} 位于疾病关键区域 "
            f"（{', '.join(r['name'] for r in result['matched_regions'])}），"
            f"但附近无已标注关键残基，建议持续监测或功能验证。"
        )
    elif nearby:
        result["upgrade_recommendation"] = "monitor"
        result["reasoning"] = (
            f"{gene} p.{protein_position} 邻近关键残基 "
            f"（{', '.join(f'{n['residue']}({n['distance']}aa)' for n in nearby[:3])}），"
            f"但不在已标注的关键区域内。"
        )
    else:
        result["reasoning"] = (
            f"{gene} p.{protein_position} 不在 {disease_name} 已标注的 "
            f"疾病相关区域或关键残基附近。"
        )

    return result


def run_domain_dive_for_variants(
    variants: list[dict],
    disease_name: str,
    min_tier: int = 2,
    max_tier: int = 3,
) -> list[dict]:
    """Run domain-dive on a list of variant dicts.

    Only variants whose 'tier' is within [min_tier, max_tier] are analyzed.
    No AF or ClinVar pre-filtering is applied: domain-dive is a targeted
    structural/contextual exploration for core genes defined in the disease
    template. Each variant dict is expected to have 'gene', 'hgvsp', and
    optionally 'tier'.
    """
    results: list[dict] = []
    for v in variants:
        tier = v.get("tier", 3)
        if tier < min_tier or tier > max_tier:
            continue
        gene = v.get("gene")
        hgvsp = v.get("hgvsp")
        pos = parse_protein_position(hgvsp)
        if not gene or pos is None:
            continue
        assessment = assess_variant_in_key_regions(gene, pos, disease_name, hgvsp)
        # Only return variants with some evidence
        if assessment["upgrade_recommendation"] != "no_evidence":
            results.append({"variant": v, "domain_dive": assessment})
    return results


if __name__ == "__main__":
    # Quick self-test with XDH p.Asn1087Ser
    import json
    test = assess_variant_in_key_regions("XDH", 1087, "高尿酸", "p.Asn1087Ser")
    print(json.dumps(test, indent=2, ensure_ascii=False))
