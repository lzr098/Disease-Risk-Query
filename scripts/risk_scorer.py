"""Composite risk scoring from GPA tiers, literature support, and rarity."""

from __future__ import annotations

import logging
from typing import Any, Optional

from constants import (
    COMPLEX_GWAS_HET_POINTS,
    COMPLEX_GWAS_HOM_POINTS,
    COMPLEX_GWAS_MAX_POINTS,
    COMPLEX_LIT_GENE_BONUS,
    COMPLEX_LIT_VARIANT_BONUS,
    COMPLEX_RARE_CLINVAR_BENIGN_PENALTY,
    COMPLEX_RARE_CLINVAR_LP_P_POINTS,
    COMPLEX_RARE_LOF_BASE,
    COMPLEX_RARE_MAX_PER_VARIANT,
    COMPLEX_RARE_MISSENSE_BASE,
    COMPLEX_WEIGHT_GWAS_COMMON,
    COMPLEX_WEIGHT_LITERATURE,
    COMPLEX_WEIGHT_MONOGENIC,
    COMPLEX_WEIGHT_RARE_FUNCTIONAL,
    COMPLEX_WEIGHT_RARITY,
    CONTRIBUTION_HIGH_THRESHOLD,
    CONTRIBUTION_LOW_THRESHOLD,
    CONTRIBUTION_MODERATE_THRESHOLD,
    LIT_GENE_BONUS,
    LIT_VARIANT_BONUS,
    RISK_HIGH_THRESHOLD,
    RISK_LOW_THRESHOLD,
    RISK_MODERATE_THRESHOLD,
    TIER1_BASE,
    TIER1_HOM_PURE_LOF_BONUS,
    TIER2_BASE,
    TIER2_MULTIHIT_BONUS,
    VCF_GWAS_DOWNWEIGHT_IF_FILTERED,
    WEIGHT_LITERATURE,
    WEIGHT_RARITY,
    WEIGHT_SEX_AGE,
    WEIGHT_TIER1,
    WEIGHT_TIER2,
)

logger = logging.getLogger(__name__)

HIGH_IMPACT = {"HIGH"}
LOF_CONSEQUENCES = {
    "stop_gained",
    "frameshift_variant",
    "splice_donor_variant",
    "splice_acceptor_variant",
    "start_lost",
    "stop_lost",
}


def _is_pure_lof(variant: dict) -> bool:
    impact = (variant.get("IMPACT") or "").upper()
    cons = (variant.get("Consequence") or "").lower()
    return impact in HIGH_IMPACT and any(c in cons for c in LOF_CONSEQUENCES)


def _is_homozygous(variant: dict) -> bool:
    gt = variant.get("GT") or variant.get("genotype") or ""
    return gt in {"1/1", "1|1"}


def _gnomad_af(variant: dict) -> Optional[float]:
    for key in ("gnomad_af", "gnomAD_AF", "AF", "af"):
        val = variant.get(key)
        if val is not None and val != "":
            try:
                return float(val)
            except (TypeError, ValueError):
                continue
    return None


def _variant_id(variant: dict) -> str:
    """Build chrom:pos:ref:alt from GPA lowercase or VCF uppercase keys."""
    chrom = variant.get("chrom") or variant.get("CHROM") or ""
    pos = variant.get("pos") or variant.get("POS") or ""
    ref = variant.get("ref") or variant.get("REF") or ""
    alt = variant.get("alt") or variant.get("ALT") or ""
    return f"{chrom}:{pos}:{ref}:{alt}"


def _variant_gene(variant: dict) -> str:
    return variant.get("gene") or variant.get("GENE") or ""


def score_variants(
    tier1: list[dict],
    tier2: list[dict],
    multi_hit_genes: list[str],
    literature_variants: Optional[list[dict]] = None,
    literature_genes: Optional[set[str]] = None,
    gene_contribution_map: Optional[dict[str, float]] = None,
) -> dict:
    """Score individual variant contributions.

    If gene_contribution_map is provided, per-variant base contribution is
    multiplied by the gene's contribution_score (capped at 1.0).

    Returns dict with per-variant scores and component sums.
    """
    gene_contribution_map = gene_contribution_map or {}
    literature_variants = literature_variants or []
    literature_genes = literature_genes or set()

    lit_variant_set: set[str] = set()
    for lv in literature_variants:
        chrom = lv.get("CHROM") or lv.get("chrom")
        pos = lv.get("POS") or lv.get("pos")
        ref = lv.get("REF") or lv.get("ref")
        alt = lv.get("ALT") or lv.get("alt")
        if chrom and pos:
            lit_variant_set.add(f"{chrom}:{pos}:{ref}:{alt}")

    tier1_score = 0
    tier1_details = []
    for v in tier1:
        gene = _variant_gene(v)
        gene_weight = min(1.0, gene_contribution_map.get(gene, 0.5) + 0.5)
        contribution = round(TIER1_BASE * gene_weight, 1)
        flags = []
        if _is_pure_lof(v) and _is_homozygous(v):
            contribution += TIER1_HOM_PURE_LOF_BONUS
            flags.append("homozygous_LoF")
        tier1_score += contribution
        v_key = _variant_id(v)
        if v_key in lit_variant_set:
            flags.append("LITERATURE_SUPPORTED")
        tier1_details.append({
            "gene": gene,
            "variant": v_key,
            "contribution": contribution,
            "gene_contribution_score": gene_contribution_map.get(gene),
            "flags": flags,
        })

    tier2_score = 0
    tier2_details = []
    for v in tier2:
        gene = _variant_gene(v)
        gene_weight = min(1.0, gene_contribution_map.get(gene, 0.3) + 0.5)
        contribution = round(TIER2_BASE * gene_weight, 1)
        flags = []
        if gene in multi_hit_genes:
            contribution += TIER2_MULTIHIT_BONUS
            flags.append("multi_hit")
        tier2_score += contribution
        v_key = _variant_id(v)
        if v_key in lit_variant_set:
            flags.append("LITERATURE_SUPPORTED")
        elif gene in literature_genes:
            flags.append("LITERATURE_GENE_SUPPORTED")
        tier2_details.append({
            "gene": gene,
            "variant": v_key,
            "contribution": contribution,
            "gene_contribution_score": gene_contribution_map.get(gene),
            "flags": flags,
        })

    # Literature component: count distinct variant-level and gene-level supports
    lit_variant_hits = sum(
        1 for d in tier1_details + tier2_details if "LITERATURE_SUPPORTED" in d["flags"]
    )
    lit_gene_hits = sum(
        1 for d in tier2_details if "LITERATURE_GENE_SUPPORTED" in d["flags"]
    )
    lit_score = min(WEIGHT_LITERATURE, lit_variant_hits * LIT_VARIANT_BONUS + lit_gene_hits * LIT_GENE_BONUS)

    # Rarity component: use lowest AF among tier1/2 variants
    all_scored = tier1 + tier2
    afs = [_gnomad_af(v) for v in all_scored if _gnomad_af(v) is not None]
    min_af = min(afs) if afs else None
    if min_af is None:
        rarity_score = WEIGHT_RARITY  # conservative: assume rare
    elif min_af < 0.001:
        rarity_score = WEIGHT_RARITY
    elif min_af < 0.01:
        rarity_score = WEIGHT_RARITY * 0.7
    elif min_af < 0.05:
        rarity_score = WEIGHT_RARITY * 0.3
    else:
        rarity_score = 0

    return {
        "tier1_score": tier1_score,
        "tier1_details": tier1_details,
        "tier2_score": tier2_score,
        "tier2_details": tier2_details,
        "literature_score": lit_score,
        "rarity_score": rarity_score,
        "min_gnomad_af": min_af,
    }


def calculate_total_score(
    tier1: list[dict],
    tier2: list[dict],
    multi_hit_genes: list[str],
    literature_variants: Optional[list[dict]] = None,
    literature_genes: Optional[set[str]] = None,
    sex_age_bonus: int = 0,
    gene_contribution_map: Optional[dict[str, float]] = None,
) -> dict:
    """Calculate weighted total score (0-100)."""
    components = score_variants(
        tier1, tier2, multi_hit_genes,
        literature_variants=literature_variants,
        literature_genes=literature_genes,
        gene_contribution_map=gene_contribution_map,
    )

    raw_total = (
        min(WEIGHT_TIER1, components["tier1_score"]) +
        min(WEIGHT_TIER2, components["tier2_score"]) +
        components["literature_score"] +
        components["rarity_score"] +
        min(WEIGHT_SEX_AGE, sex_age_bonus)
    )
    total = min(100, max(0, raw_total))

    if total >= RISK_HIGH_THRESHOLD:
        level = "高风险"
        meaning = "携带明确致病突变，强烈建议遗传咨询"
    elif total >= RISK_MODERATE_THRESHOLD:
        level = "中等风险"
        meaning = "携带可能致病或多证据支持变异，建议监测"
    elif total >= RISK_LOW_THRESHOLD:
        level = "低风险"
        meaning = "有变异但证据较弱，或仅为携带者"
    else:
        level = "无明确风险"
        meaning = "未检出相关变异"

    return {
        "total_score": round(total, 1),
        "risk_level": level,
        "risk_meaning": meaning,
        "components": components,
        "sex_age_bonus": sex_age_bonus,
    }


def sex_age_bonus(
    disease_name: str,
    sex: str,
    age: Optional[int],
    genes: list[str],
) -> int:
    """Compute sex/age correction bonus.

    Simple rule: if disease name or gene list hints at X-linked inheritance
    and sex matches male for X-linked recessive, add small bonus. Age below
    typical onset lowers risk by negative value (caller should interpret).
    """
    low = disease_name.lower()
    bonus = 0
    x_linked_keywords = {"x-linked", "xlinked", "hemophilia", "duchenne", "fragile x"}
    if any(kw in low for kw in x_linked_keywords):
        if sex.lower() == "male":
            bonus += 5
    return min(bonus, 5)


# ---------------------------------------------------------------------------
# Complex-trait contribution scoring
# ---------------------------------------------------------------------------


def _is_likely_benign_clinvar(variant: dict) -> bool:
    sig = (variant.get("CLIN_SIG") or variant.get("clinvar") or "").lower()
    return "benign" in sig and "pathogenic" not in sig


def _is_likely_pathogenic_clinvar(variant: dict) -> bool:
    sig = (variant.get("CLIN_SIG") or variant.get("clinvar") or "").lower()
    return "pathogenic" in sig and "benign" not in sig


def _consequence_contains(variant: dict, term: str) -> bool:
    cons = (variant.get("Consequence") or variant.get("consequence") or "").lower()
    return term in cons


def _is_missense(variant: dict) -> bool:
    return _consequence_contains(variant, "missense")


def score_rare_functional_variants(
    variants: list[dict],
    core_genes: set[str],
    literature_variants: Optional[list[dict]] = None,
    literature_genes: Optional[set[str]] = None,
    gene_contribution_map: Optional[dict[str, float]] = None,
) -> dict:
    """Score rare functional variants for complex-trait contribution.

    Each variant gets a contribution score (0-20) based on:
      - gene membership in core set
      - predicted functional impact (LoF vs missense)
      - ClinVar interpretation
      - zygosity
      - population rarity (gnomAD AF)

    Returns a dict with per-variant details and the capped dimension score.
    """
    literature_variants = literature_variants or []
    literature_genes = literature_genes or set()
    gene_contribution_map = gene_contribution_map or {}

    lit_variant_set: set[str] = set()
    for lv in literature_variants:
        chrom = lv.get("CHROM") or lv.get("chrom")
        pos = lv.get("POS") or lv.get("pos")
        ref = lv.get("REF") or lv.get("ref")
        alt = lv.get("ALT") or lv.get("alt")
        if chrom and pos:
            lit_variant_set.add(f"{chrom}:{pos}:{ref}:{alt}")

    details = []
    total = 0.0
    for v in variants:
        gene = _variant_gene(v)
        if gene.upper() not in {g.upper() for g in core_genes}:
            continue
        if _is_likely_benign_clinvar(v):
            continue

        points = 0.0
        notes = []
        if _is_pure_lof(v):
            points += COMPLEX_RARE_LOF_BASE
            notes.append("LoF")
        elif _is_missense(v):
            points += COMPLEX_RARE_MISSENSE_BASE
            notes.append("missense")
        else:
            # Other moderate/high-impact variants get a small base score
            impact = (v.get("IMPACT") or v.get("impact") or "").upper()
            if impact in {"HIGH", "MODERATE"}:
                points += 2
                notes.append(f"{impact.lower()}_impact")
            else:
                continue

        if _is_likely_pathogenic_clinvar(v):
            points += COMPLEX_RARE_CLINVAR_LP_P_POINTS
            notes.append("ClinVar_LP/P")

        af = _gnomad_af(v)
        if af is not None:
            if af < 0.0001:
                points += 4
                notes.append("AF<0.01%")
            elif af < 0.001:
                points += 3
                notes.append("AF<0.1%")
            elif af < 0.01:
                points += 1
                notes.append("AF<1%")

        if _is_homozygous(v):
            points += 3
            notes.append("homozygous")

        v_key = _variant_id(v)
        if v_key in lit_variant_set:
            points += COMPLEX_LIT_VARIANT_BONUS
            notes.append("literature_variant")
        elif gene.upper() in {g.upper() for g in literature_genes}:
            points += COMPLEX_LIT_GENE_BONUS
            notes.append("literature_gene")

        # Apply gene contribution weight from disease template
        gene_weight = min(1.0, gene_contribution_map.get(gene, 0.2) + 0.5)
        points = points * gene_weight

        points = min(COMPLEX_RARE_MAX_PER_VARIANT, max(0, points))
        total += points
        details.append({
            "gene": gene,
            "variant": v_key,
            "contribution": round(points, 1),
            "gene_contribution_score": gene_contribution_map.get(gene),
            "notes": notes,
            "af": af,
        })

    dimension_score = min(COMPLEX_WEIGHT_RARE_FUNCTIONAL, total)
    return {
        "score": round(dimension_score, 1),
        "raw_score": round(total, 1),
        "details": details,
    }


def score_gwas_contributions(
    gwas_lead_snps: list[dict],
    snp_contribution_map: Optional[dict[str, float]] = None,
) -> dict:
    """Score direct GWAS lead SNP genotypes for complex traits.

    If snp_contribution_map is provided, each hit is weighted by the SNP's
    contribution_score (or abs(beta) when available) instead of flat points.
    Heterozygous risk allele uses full weight; homozygous doubles it.
    Capped at COMPLEX_GWAS_MAX_POINTS.
    """
    snp_contribution_map = snp_contribution_map or {}
    total = 0.0
    hits = []
    for s in gwas_lead_snps:
        gt = s.get("sample_gt") or {}
        genotype = gt.get("gt", "")
        if not genotype:
            continue
        rsid = s.get("rsid", "")
        # Dosage
        if genotype in {"1/1", "1|1"}:
            dosage = 2
        elif genotype in {"0/1", "0|1", "1/0", "1|0"}:
            dosage = 1
        else:
            continue

        # Weighted contribution from template
        weight = snp_contribution_map.get(rsid)
        if weight is None:
            beta = s.get("beta")
            if beta is not None:
                weight = min(0.30, max(0.03, abs(beta) * 1.2))
            else:
                weight = 0.10
        points = dosage * weight * 10  # scale to roughly 0-6 per hit
        total += points
        hits.append({
            "rsid": rsid,
            "gene": s.get("gene"),
            "chrom": gt.get("chrom"),
            "pos": gt.get("pos"),
            "gt": genotype,
            "dosage": dosage,
            "snp_weight": round(weight, 3),
            "points": round(points, 2),
        })

    dimension_score = min(COMPLEX_GWAS_MAX_POINTS, total)
    return {
        "score": round(dimension_score, 1),
        "raw_score": round(total, 1),
        "hits": hits,
    }


def score_literature_contributions(
    tier1: list[dict],
    tier2: list[dict],
    tier3: list[dict],
    literature_variants: Optional[list[dict]] = None,
    literature_genes: Optional[set[str]] = None,
) -> dict:
    """Score literature support across all tiers for complex traits."""
    literature_variants = literature_variants or []
    literature_genes = literature_genes or set()

    lit_variant_set: set[str] = set()
    for lv in literature_variants:
        chrom = lv.get("CHROM") or lv.get("chrom")
        pos = lv.get("POS") or lv.get("pos")
        ref = lv.get("REF") or lv.get("ref")
        alt = lv.get("ALT") or lv.get("alt")
        if chrom and pos:
            lit_variant_set.add(f"{chrom}:{pos}:{ref}:{alt}")

    total = 0.0
    details = []
    for v in tier1 + tier2 + tier3:
        gene = _variant_gene(v)
        v_key = _variant_id(v)
        if v_key in lit_variant_set:
            total += COMPLEX_LIT_VARIANT_BONUS
            details.append({"gene": gene, "variant": v_key, "type": "variant", "points": COMPLEX_LIT_VARIANT_BONUS})
        elif gene.upper() in {g.upper() for g in literature_genes}:
            total += COMPLEX_LIT_GENE_BONUS
            details.append({"gene": gene, "variant": v_key, "type": "gene", "points": COMPLEX_LIT_GENE_BONUS})

    dimension_score = min(COMPLEX_WEIGHT_LITERATURE, total)
    return {
        "score": round(dimension_score, 1),
        "raw_score": round(total, 1),
        "details": details,
    }


def score_rarity_contribution(
    variants: list[dict],
) -> dict:
    """Score aggregate rarity/protective evidence for complex traits.

    Uses the lowest gnomAD AF among scored variants; rare variants support a
    stronger genetic contribution. Protective allele logic can be added later.
    """
    afs = [_gnomad_af(v) for v in variants if _gnomad_af(v) is not None]
    min_af = min(afs) if afs else None
    if min_af is None:
        score = COMPLEX_WEIGHT_RARITY
        note = "no AF data; conservatively assume rare"
    elif min_af < 0.0001:
        score = COMPLEX_WEIGHT_RARITY
        note = "min AF < 0.01%"
    elif min_af < 0.001:
        score = COMPLEX_WEIGHT_RARITY * 0.7
        note = "min AF < 0.1%"
    elif min_af < 0.01:
        score = COMPLEX_WEIGHT_RARITY * 0.3
        note = "min AF < 1%"
    else:
        score = 0.0
        note = "min AF >= 1%"

    return {
        "score": round(score, 1),
        "min_gnomad_af": min_af,
        "note": note,
    }


def score_monogenic_contributions(
    tier1: list[dict],
    tier2: list[dict],
) -> dict:
    """Score high-penetrance monogenic contributions in complex traits.

    Tier 1 variants are counted at full weight, Tier 2 at half weight.
    Capped at COMPLEX_WEIGHT_MONOGENIC.
    """
    total = 0.0
    details = []
    for v in tier1:
        total += COMPLEX_WEIGHT_MONOGENIC
        details.append({
            "gene": _variant_gene(v),
            "variant": _variant_id(v),
            "tier": 1,
            "points": COMPLEX_WEIGHT_MONOGENIC,
        })
    for v in tier2:
        total += COMPLEX_WEIGHT_MONOGENIC / 2
        details.append({
            "gene": _variant_gene(v),
            "variant": _variant_id(v),
            "tier": 2,
            "points": COMPLEX_WEIGHT_MONOGENIC / 2,
        })

    dimension_score = min(COMPLEX_WEIGHT_MONOGENIC, total)
    return {
        "score": round(dimension_score, 1),
        "raw_score": round(total, 1),
        "details": details,
    }


def calculate_contribution_score(
    tier1: list[dict],
    tier2: list[dict],
    tier3: list[dict],
    core_genes: set[str],
    gwas_lead_snps: Optional[list[dict]] = None,
    literature_variants: Optional[list[dict]] = None,
    literature_genes: Optional[set[str]] = None,
    vcf_qc: Optional[dict] = None,
    gene_contribution_map: Optional[dict[str, float]] = None,
    snp_contribution_map: Optional[dict[str, float]] = None,
) -> dict:
    """Calculate complex-trait genetic contribution score (0-100).

    Returns a structured result with per-dimension and per-variant contributions,
    suitable for reporting "how much might genetics contribute to this phenotype".

    If the VCF appears to have been hard-filtered for common variants (low
    presence rate of GWAS anchor SNPs), the GWAS dimension is downweighted and
    the total is rescaled so that scores remain comparable across VCF types.
    """
    gwas_lead_snps = gwas_lead_snps or []
    vcf_qc = vcf_qc or {}

    monogenic = score_monogenic_contributions(tier1, tier2)
    rare = score_rare_functional_variants(
        tier3, core_genes,
        literature_variants=literature_variants,
        literature_genes=literature_genes,
        gene_contribution_map=gene_contribution_map,
    )
    gwas = score_gwas_contributions(gwas_lead_snps, snp_contribution_map=snp_contribution_map)
    lit = score_literature_contributions(
        tier1, tier2, tier3,
        literature_variants=literature_variants,
        literature_genes=literature_genes,
    )
    rarity = score_rarity_contribution(tier1 + tier2 + tier3)

    # Adjust GWAS dimension if VCF appears to have filtered common variants
    gwas_adjustment = None
    if vcf_qc.get("common_variants_filtered"):
        effective_gwas_weight = COMPLEX_WEIGHT_GWAS_COMMON * VCF_GWAS_DOWNWEIGHT_IF_FILTERED
        original_gwas_score = gwas["score"]
        gwas["original_score"] = original_gwas_score
        gwas["score"] = round(min(effective_gwas_weight, original_gwas_score), 1)
        unassessed_gwas = round(COMPLEX_WEIGHT_GWAS_COMMON - effective_gwas_weight, 1)
        gwas_adjustment = {
            "downweight_factor": VCF_GWAS_DOWNWEIGHT_IF_FILTERED,
            "effective_gwas_weight": round(effective_gwas_weight, 1),
            "unassessed_gwas_weight": unassessed_gwas,
            "presence_rate": vcf_qc.get("presence_rate", 0.0),
            "note": (
                f"VCF 疑似过滤常见变异（锚定位点检出率 {vcf_qc.get('presence_rate', 0):.0%}），"
                f"GWAS 维度按 {VCF_GWAS_DOWNWEIGHT_IF_FILTERED:.0%} 折算，"
                f"{unassessed_gwas} 分标记为未评估。"
            ),
        }
        gwas["adjustment"] = gwas_adjustment

    unassessed = gwas_adjustment["unassessed_gwas_weight"] if gwas_adjustment else 0.0
    effective_max = 100.0 - unassessed

    raw_total = (
        monogenic["score"]
        + rare["score"]
        + gwas["score"]
        + lit["score"]
        + rarity["score"]
    )
    # Rescale to 0-100 so scores remain comparable across filtered/unfiltered VCFs
    total = raw_total * (100.0 / effective_max) if effective_max > 0 else raw_total
    total = min(100.0, max(0.0, total))

    if total >= CONTRIBUTION_HIGH_THRESHOLD:
        level = "遗传贡献较高"
        meaning = (
            "遗传因素在当前表型中占比较重；建议结合家族史和长期生化指标关注。"
        )
    elif total >= CONTRIBUTION_MODERATE_THRESHOLD:
        level = "遗传贡献中等"
        meaning = "遗传可能是多个致病因素之一，生活方式和代谢状态仍占重要比例。"
    elif total >= CONTRIBUTION_LOW_THRESHOLD:
        level = "遗传贡献较低"
        meaning = "检出少量相关变异，但不足以单独解释当前表型。"
    else:
        level = "无明确遗传贡献"
        meaning = "未检出对该表型有实质影响的遗传证据。"

    result = {
        "mode": "complex",
        "total_score": round(total, 1),
        "contribution_level": level,
        "contribution_meaning": meaning,
        "components": {
            "monogenic_score": monogenic["score"],
            "monogenic_details": monogenic["details"],
            "rare_functional_score": rare["score"],
            "rare_functional_details": rare["details"],
            "gwas_score": gwas["score"],
            "gwas_hits": gwas["hits"],
            "literature_score": lit["score"],
            "literature_details": lit["details"],
            "rarity_score": rarity["score"],
            "rarity_note": rarity["note"],
            "min_gnomad_af": rarity["min_gnomad_af"],
        },
    }
    if gwas_adjustment:
        result["gwas_adjustment"] = gwas_adjustment
        result["unassessed_weight"] = unassessed
        result["raw_total_before_rescale"] = round(raw_total, 1)
        result["components"]["gwas_original_score"] = gwas.get("original_score", gwas["score"])
    return result
