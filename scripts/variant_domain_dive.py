"""Targeted variant domain-dive analysis.

Given a Tier 3 (or borderline) variant, check whether it falls into a
 disease-specific key functional region or near a critical residue.
This provides structured evidence for manual Tier 2 upgrade decisions.

v2.1: Added protein context assessment — auto-detects disordered regions,
phosphorylation sites, and known domains from VEP/UniProt feature data
to provide structural context for contribution adjustments.
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


def assess_protein_context(
    protein_position: int,
    vep_features: Optional[dict] = None,
    uniprot_features: Optional[list[dict]] = None,
) -> dict[str, Any]:
    """Assess protein structural context of a missense variant position.

    Uses VEP/UniProt feature annotations to determine whether a position
    falls in a disordered region, known domain, or phosphorylation site.

    Args:
        protein_position: 1-based amino acid position
        vep_features: VEP DOMAINS annotation (parsed, list of dicts with
                      'db', 'name', 'start', 'end')
        uniprot_features: UniProt feature list from API (list of dicts with
                         'type', 'description', 'begin', 'end')

    Returns dict with:
      - is_disordered: bool
      - is_phosphosite: bool
      - is_predictive_phosphosite: bool (adjacent Ser/Thr/Tyr within ±2)
      - domain_names: list[str] — names of domains covering this position
      - matched_features: list[str] — human-readable feature descriptions
      - downgrade_factor: float — suggested contribution multiplier
        (1.0=no downgrade, 0.5=disordered+no PTM, 0.7=disordered phosphosite)
      - reasoning: str — human-readable explanation
    """
    result: dict[str, Any] = {
        "protein_position": protein_position,
        "is_disordered": False,
        "is_phosphosite": False,
        "is_predictive_phosphosite": False,
        "domain_names": [],
        "matched_features": [],
        "downgrade_factor": 1.0,
        "reasoning": "",
    }

    disorder_keywords = [
        "disordered", "disorder", "unstructured",
        "intrinsically_disordered", "flexible_loop",
        "low_complexity", "compositional_bias",
    ]
    phospho_keywords = ["phospho", "phosphoryl"]
    domain_exclude = [
        "signal_peptide", "signal", "transit_peptide",
        "propeptide", "chain", "mature_chain",
    ]

    all_features: list[dict] = []
    if vep_features:
        for feat in vep_features:
            feat["source"] = "vep"
            all_features.append(feat)
    if uniprot_features:
        for feat in uniprot_features:
            feat["source"] = "uniprot"
            all_features.append(feat)

    for feat in all_features:
        start = feat.get("start", feat.get("begin", 0))
        end = feat.get("end", finish if (finish := feat.get("finish")) else 0)
        if not (start and end):
            continue
        if not (start <= protein_position <= end):
            continue

        description = (feat.get("description") or feat.get("name") or "").lower()
        feat_type = (feat.get("type") or "").lower()

        # Check disordered
        is_disorder = any(kw in description or kw in feat_type for kw in disorder_keywords)
        if is_disorder:
            result["is_disordered"] = True
            result["matched_features"].append(
                f"Disordered region [{start}-{end}]: {feat.get('description', feat.get('name', ''))}"
            )

        # Check phosphorylation
        is_phospho = any(kw in description or kw in feat_type for kw in phospho_keywords)
        if is_phospho:
            result["is_phosphosite"] = True
            result["matched_features"].append(
                f"Phosphosite [{start}-{end}]: {feat.get('description', feat.get('name', ''))}"
            )

        # Collect domain names
        feat_name = feat.get("name", feat.get("description", ""))
        is_domain = (
            feat.get("db") not in (None, "")
            or feat.get("source") == "vep"
            and feat_name
        )
        if is_domain and feat_name:
            name_lower = feat_name.lower()
            if not any(excl in name_lower for excl in domain_exclude):
                already = any(d.lower() == name_lower for d in result["domain_names"])
                if not already:
                    result["domain_names"].append(feat_name)

        # Check if not in any functional domain
        if feat_type in ("domain", "region", "repeat", "motif",
                         "coiled_coil", "zinc_finger", "helix", "strand",
                         "topological_domain", "transmembrane"):
            feat_name_val = feat.get("name", feat.get("description", ""))
            name_lower = feat_name_val.lower()
            if not any(excl in name_lower for excl in domain_exclude):
                already = any(d.lower() == feat_name_val.lower() for d in result["domain_names"])
                if not already:
                    result["domain_names"].append(feat_name_val)

    # Predictive phosphosite: nearby Ser/Thr/Tyr (+-2)
    # This is only checked if VEP-provided protein features are scarce
    if not result["is_phosphosite"] and vep_features:
        aa_around = _guess_nearby_phosphosite(vep_features, protein_position)
        if aa_around:
            result["is_predictive_phosphosite"] = aa_around

    # Compute downgrade factor
    reasons = []
    if result["is_disordered"]:
        if result["is_phosphosite"]:
            result["downgrade_factor"] = 0.7
            reasons.append("无序区磷酸化位点 — 功能重要性较高，部分降权(×0.7)")
        else:
            result["downgrade_factor"] = 0.5
            reasons.append("位于无序区且无PTM — 功能性证据弱，显著降权(×0.5)")
    elif not result["domain_names"] and vep_features and len(vep_features) > 0:
        # In protein but not in any known domain
        result["downgrade_factor"] = 0.8
        reasons.append("蛋白内无已知功能域覆盖 — 轻度降权(×0.8)")

    result["reasoning"] = "; ".join(reasons) if reasons else "位于已知蛋白功能域内"
    return result


def _guess_nearby_phosphosite(
    vep_features: list[dict], position: int
) -> bool:
    """Check if the variant position is Ser/Thr/Tyr but not annotated as phosphosite.

    This is a rough heuristic: if there's no phosphosite annotation at this
    position but the surrounding protein features suggest regulatory potential,
    return True as a predictive flag.
    """
    # This requires AA context from HGVSp which we don't always have
    return False


if __name__ == "__main__":
    # Quick self-test with XDH p.Asn1087Ser
    import json
    test = assess_variant_in_key_regions("XDH", 1087, "高尿酸", "p.Asn1087Ser")
    print(json.dumps(test, indent=2, ensure_ascii=False))
