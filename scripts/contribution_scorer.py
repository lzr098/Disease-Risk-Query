"""Disease contribution scoring model.

Scores variants by how they contribute to the target disease, not just by
generic pathogenicity. Supports multiple layers:
- mendelian_high: high-penetrance pathogenic variants (from disease-space VCF)
- mendelian_mod: moderate-penetrance variants in disease genes
- known_pathogenic: curated known pathogenic variants from disease template
  (queried exactly, scored by zygosity × penetrance-adjusted contribution)
- dosage_risk: risk alleles with copy-number effect (e.g. APOE e4)
- gwas_prs: weighted polygenic contribution from common variants
  (normalised by sqrt(variant_count) for cross-disease comparability)
- regulatory: weak evidence from regulatory regions

ClinGen validity adjustment: Mendelian gene weights are modulated by
clingen_validity (Definitive ×1.0, Strong ×0.85, Moderate ×0.70,
Limited ×0.50, Disputed/Refuted ×0.20, missing ×0.65).

Key domain bonus: variants in genes with known key_domains receive a
×1.15 contribution multiplier when the variant is a coding change
(missense/inframe), reflecting increased likelihood of functional impact.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field
from typing import Any, Optional

from disease_profile import DiseaseProfile, VariantWeight
from disease_space_query import KnownVariantGenotype
from constants import VCF_GWAS_DOWNWEIGHT_IF_FILTERED

logger = logging.getLogger(__name__)


def _variant_gene(v: dict) -> str:
    return v.get("gene") or v.get("GENE") or v.get("symbol") or ""


def _variant_key(v: dict) -> str:
    """Return a chrom:pos:ref>alt string from a GPA-style variant dict."""
    chrom = v.get("chrom") if v.get("chrom") is not None else v.get("CHROM")
    pos = v.get("pos") if v.get("pos") is not None else v.get("POS")
    ref = v.get("ref") if v.get("ref") is not None else v.get("REF")
    alt = v.get("alt") if v.get("alt") is not None else v.get("ALT")
    if chrom is not None and pos is not None and ref is not None and alt is not None:
        return f"{chrom}:{pos}:{ref}>{alt}"
    return v.get("variant_id") or ""


def _gt(kg: KnownVariantGenotype) -> str:
    return kg.gt if kg.gt else "0/0"


def _gnomad_af(v: dict) -> Optional[float]:
    for key in ("gnomad_af", "gnomad_eas_af", "gnomad_nfe_af", "AF", "EAS_AF"):
        val = v.get(key)
        if isinstance(val, (int, float)) and val >= 0:
            return float(val)
    return None


def _is_high_impact(v: dict) -> bool:
    impact = v.get("impact", "")
    if isinstance(impact, str):
        return impact.upper() in ("HIGH", "MODERATE")
    return False


def _is_confident_pathogenic(v: dict) -> bool:
    """Return True if ClinVar classifies this as pathogenic or likely pathogenic."""
    sig = (v.get("clinvar_sig") or "").lower()
    if not sig or sig in ("not_provided", "not_specified", "unknown"):
        return False
    if "not" in sig:
        return False
    return "pathogenic" in sig


def _is_vus_clinvar(v: dict) -> bool:
    """Return True if ClinVar classification is VUS (uncertain significance)."""
    sig = (v.get("clinvar_sig") or "").lower()
    if "uncertain" in sig:
        return True
    # Conflicting classifications without a clear lean
    if "conflicting" in sig:
        # conflicting + benign/likely_benign → treat as VUS-equivalent
        return True
    return False


def _is_likely_benign_clinvar(v: dict) -> bool:
    sig = (v.get("clinvar_sig") or "").lower()
    return "benign" in sig and "pathogenic" not in sig


def _clinvar_vus_factor(v: dict) -> float:
    """Return contribution multiplier for ClinVar VUS.
    
    VUS = 0.5 (uncertain), Conflicting-VUS = 0.7 (some submissions disagree),
    Conflicting with benign lean = 0.4 (some evidence against pathogenicity).
    Non-VUS = 1.0 (no penalty).
    """
    if not _is_vus_clinvar(v):
        return 1.0
    sig = (v.get("clinvar_sig") or "").lower()
    # Conflicting WITH benign in the mix → lean toward reduced weight
    if "conflicting" in sig:
        if "benign" in sig and "pathogenic" not in sig.replace("conflicting_classifications_of_pathogenicity", "").replace("conflicting_interpretations_of_pathogenicity", ""):
            return 0.4
        return 0.7
    # Standard VUS
    return 0.5


def _protein_context_factor(v: dict) -> tuple[float, str]:
    """Assess protein structural context from GPA domain_info and VEP/UniProt.

    Returns (multiplier, reason_string).
    
    Uses GPA's 'domain_info' field (pre-computed by GPA core) which provides
    structured domain annotations including:
    - 'domain': domain name (e.g., 'MutS_V', 'Disordered', 'Kinase')
    - 'domain_range': residue range
    - 'function': domain function category
    
    Also checks VEP DOMAINS and UniProt features when available.
    
    Disordered region + no PTM → ×0.5
    Disordered region + phosphosite → ×0.7
    No known domain ("inter-domain") → ×0.8
    Otherwise → ×1.0
    """
    # Check GPA's domain_info first (most reliable source)
    domain_info = v.get("domain_info")
    if isinstance(domain_info, dict):
        domain = (domain_info.get("domain") or "").lower()
        function_type = (domain_info.get("function") or "").lower()
        
        if "disordered" in domain or "disordered" in function_type:
            # Check if this is also a known modification site
            note = domain_info.get("note", "").lower()
            if any(kw in note for kw in ["phospho", "modification", "acetyl"]):
                return 0.7, f"disordered region with modification site (×0.7): {domain_info.get('domain_range','?')}"
            return 0.5, f"disordered region (×0.5): {domain_info.get('domain_range','?')}"
        
        if "inter-domain" in domain or "unannotated" in domain:
            return 0.8, "inter-domain/unannotated position (×0.8)"
    
    # Fallback to VEP/UniProt features
    domains_raw = v.get("domains")
    uniprot_features = v.get("uniprot_features")
    
    if not domains_raw and not uniprot_features:
        return 1.0, ""
    
    # Extract protein position from HGVSp
    hgvsp = v.get("hgvsp") or v.get("HGVSp") or ""
    import re
    pos_match = re.search(r"p?\.[A-Za-z*]+\d+", hgvsp)
    if not pos_match:
        return 1.0, ""
    digits = re.search(r"(\d+)", pos_match.group(0))
    if not digits:
        return 1.0, ""
    protein_position = int(digits.group(1))
    
    # VEP DOMAINS: parse into feature list
    vep_features = []
    if isinstance(domains_raw, list):
        for d in domains_raw:
            if isinstance(d, dict):
                vep_features.append({
                    "db": d.get("db", ""),
                    "name": d.get("name", d.get("desc", "")),
                    "start": d.get("start", 0),
                    "end": d.get("end", 0),
                })
    
    # UniProt features
    uni_features = []
    if isinstance(uniprot_features, list):
        uni_features = uniprot_features
    
    # Use domain_dive assessor
    from variant_domain_dive import assess_protein_context
    ctx = assess_protein_context(protein_position, vep_features, uni_features)
    
    if ctx.get("downgrade_factor", 1.0) < 1.0:
        return ctx["downgrade_factor"], ctx.get("reasoning", "")
    return 1.0, ""


def _clingen_weight(clingen_validity: str, is_mendelian: bool) -> float:
    """Return ClinGen validity weight multiplier for a gene.

    Mendelian genes are expected to have a ClinGen rating.
    Missing rating on a Mendelian gene is penalized (×0.65).
    GWAS genes without ratings are not penalized (×1.0).
    """
    if not clingen_validity or not clingen_validity.strip():
        return 0.65 if is_mendelian else 1.0
    
    validity = clingen_validity.strip().lower()
    
    if "definitive" in validity:
        return 1.00
    if "strong" in validity:
        return 0.85
    if "moderate" in validity:
        return 0.70
    if "limited" in validity:
        return 0.50
    if "disputed" in validity or "refuted" in validity:
        return 0.20
    
    # Unknown rating string — treat same as missing
    return 0.65 if is_mendelian else 1.0


def _is_coding_variant(v: dict) -> bool:
    """Check if variant is a coding change (missense, inframe, stop, etc.)."""
    consequence = (v.get("consequence") or v.get("CSQ") or "").lower()
    coding_terms = [
        "missense", "inframe", "stop_gained", "stop_lost",
        "start_lost", "frameshift", "synonymous", "coding",
        "nonsynonymous", "non_synonymous", "amino_acid",
    ]
    return any(term in consequence for term in coding_terms)


def _domain_bonus(gene: str, profile: DiseaseProfile, v: dict) -> float:
    """Return domain-hit bonus multiplier for a variant.

    If the gene has key_domains annotation and the variant is a coding
    change (missense/inframe), apply a 1.15× multiplier reflecting
    increased likelihood of hitting a critical domain.
    """
    gene_map = {g.gene: g for g in profile.gene_set}
    gw = gene_map.get(gene)
    if not gw or not gw.key_domains:
        return 1.0
    if _is_coding_variant(v):
        return 1.15
    return 1.0


@dataclass
class ContributionResult:
    mendelian_high: list[dict] = field(default_factory=list)
    mendelian_mod: list[dict] = field(default_factory=list)
    known_pathogenic: list[dict] = field(default_factory=list)
    dosage_risk: list[dict] = field(default_factory=list)
    gwas_prs: dict = field(default_factory=lambda: {
        "score": 0.0,
        "percentile": None,
        "variant_count": 0,
        "variants": [],
    })
    regulatory: list[dict] = field(default_factory=list)
    prs_high: dict = field(default_factory=lambda: {
        "score": 0.0,
        "variant_count": 0,
        "variants": [],
    })
    clinvar_enriched: list[dict] = field(default_factory=list)
    overall_level: str = "uncertain"
    overall_score: Optional[float] = None
    details: list[dict] = field(default_factory=list)
    layer_levels: dict[str, str] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "mendelian_high": self.mendelian_high,
            "mendelian_mod": self.mendelian_mod,
            "known_pathogenic": self.known_pathogenic,
            "dosage_risk": self.dosage_risk,
            "gwas_prs": self.gwas_prs,
            "regulatory": self.regulatory,
            "prs_high": self.prs_high,
            "clinvar_enriched": self.clinvar_enriched,
            "overall_level": self.overall_level,
            "overall_score": self.overall_score,
            "details": self.details,
            "layer_levels": self.layer_levels,
        }


def _score_mendelian_high(
    profile: DiseaseProfile,
    tiered: dict[str, list[dict]],
) -> list[dict]:
    """High-penetrance pathogenic variants in mendelian_high genes."""
    results = []
    gene_map = {g.gene: g for g in profile.gene_set}
    for v in tiered.get("tier1_variants", []):
        gene = _variant_gene(v)
        gw = gene_map.get(gene)
        if not gw or gw.tier != "mendelian_high":
            continue
        if _is_likely_benign_clinvar(v):
            continue
        
        clingen_mult = _clingen_weight(gw.clingen_validity, True)
        domain_mult = _domain_bonus(gene, profile, v)
        contribution = 1.0 * clingen_mult * domain_mult
        
        note_parts = ["High-penetrance variant in mendelian gene"]
        if clingen_mult < 1.0:
            note_parts.append(f"clingen_validity={gw.clingen_validity or 'missing'} (×{clingen_mult})")
        if domain_mult > 1.0:
            note_parts.append("coding change in gene with key domains (×1.15)")
        
        results.append({
            "gene": gene,
            "variant": _variant_key(v),
            "tier": 1,
            "contribution": round(min(contribution, 1.0), 3),
            "note": "; ".join(note_parts),
            "raw": v,
        })
    return results


def _is_damaging_or_pathogenic(v: dict) -> bool:
    """Return True if variant is HIGH/MODERATE impact or ClinVar P/LP."""
    if _is_high_impact(v):
        return True
    sig = (v.get("clinvar_sig") or "").lower()
    return "pathogenic" in sig


def _score_mendelian_mod(
    profile: DiseaseProfile,
    tiered: dict[str, list[dict]],
) -> list[dict]:
    """Moderate-penetrance variants in mendelian_mod or mendelian_high genes."""
    results = []
    gene_map = {g.gene: g for g in profile.gene_set}
    for tier_key in ("tier1_variants", "tier2_variants", "tier3_variants"):
        for v in tiered.get(tier_key, []):
            gene = _variant_gene(v)
            gw = gene_map.get(gene)
            if not gw or gw.tier not in ("mendelian_mod", "mendelian_high"):
                continue
            if _is_likely_benign_clinvar(v):
                continue
            # For Tier 3, require additional evidence (pathogenic or HIGH impact + rare)
            tier = 1 if tier_key == "tier1_variants" else (2 if tier_key == "tier2_variants" else 3)
            if tier == 3:
                af = _gnomad_af(v)
                if af is not None and af > 0.01:
                    continue
                if not _is_damaging_or_pathogenic(v):
                    continue
            af = _gnomad_af(v)
            rarity_bonus = 1.0 if af is None else max(0.5, 1.0 - af)
            contribution = gw.contribution_score * gw.penetrance_score * rarity_bonus
            # Tier 3 variants get a lower ceiling
            if tier == 3:
                contribution *= 0.5
            # Zygosity factor
            gt = (v.get("gt") or "./.").replace("|", "/")
            zygosity_factor = 1.5 if gt in ("1/1", "1|1") else 1.0
            contribution *= zygosity_factor
            
            # ClinGen validity and domain context
            clingen_mult = _clingen_weight(gw.clingen_validity, True)
            # VUS-like variants in core genes are domain-aware: outside known
            # functional domains they are downgraded, in-domain variants keep
            # the established key-domain bonus.
            domain_factor, domain_reason = _protein_context_factor(v)
            is_vus_like = not _is_confident_pathogenic(v)
            if is_vus_like and domain_factor < 1.0:
                contribution *= domain_factor
            # Apply key-domain bonus only for in-domain variants or confirmed P/LP
            domain_mult = 1.15 if (domain_factor >= 1.0 or not is_vus_like) and _is_coding_variant(v) else 1.0
            contribution *= clingen_mult * domain_mult
            
            note_parts = [f"{gw.penetrance} penetrance variant in {gene}"]
            if clingen_mult < 1.0:
                note_parts.append(f"clingen_validity={gw.clingen_validity or 'missing'} (×{clingen_mult})")
            if domain_factor < 1.0 and is_vus_like:
                note_parts.append(domain_reason)
            if domain_mult > 1.0:
                note_parts.append("coding change in gene with key domains (×1.15)")
            
            results.append({
                "gene": gene,
                "variant": _variant_key(v),
                "tier": tier,
                "contribution": round(min(contribution, 1.0), 3),
                "note": "; ".join(note_parts),
                "raw": v,
            })
    return results


def _score_known_pathogenic(
    profile: DiseaseProfile,
    known_genotypes: list[KnownVariantGenotype],
) -> list[dict]:
    """Curated known pathogenic variants from the disease template.

    Uses a penetrance-informed zygosity model to distinguish:
    - dominant variants (1 copy sufficient)
    - recessive variants (carrier → minimal, homozygous → full)
    - dose-dependent risk variants (homozygous amplified for moderate+ confidence)

    For low-penetrance genes (penetrance < 0.3), homozygous genotypes that are
    moderate/high confidence receive a 2.0× amplifier. This captures common
    risk variants like ABCG2 Q141K in hyperuricemia that have strong
    population-level dosage effects despite low monogenic penetrance.
    """
    results = []
    gene_map = {g.gene: g for g in profile.gene_set}
    pathogenic_vars = [
        kg for kg in known_genotypes
        if kg.variant.variant_class == "known_pathogenic"
    ]
    for kg in pathogenic_vars:
        v = kg.variant
        if kg.dosage == 0 and not kg.inferred_ref_ref:
            continue  # legitimately absent

        if kg.dosage == 0 and kg.inferred_ref_ref:
            # REF/REF inferred — record for reporting but zero contribution
            results.append({
                "rsid": v.rsid or v.vcf_key,
                "gene": v.gene or "-",
                "variant": v.vcf_key,
                "gt": kg.gt,
                "dosage": 0,
                "risk_allele": v.effect_allele or v.alt,
                "contribution": 0.0,
                "inferred_ref_ref": True,
                "confidence": v.confidence,
                "note": v.note or "",
            })
            continue
        gene_w = gene_map.get(v.gene)
        penetrance_score = gene_w.penetrance_score if gene_w else 0.5
        confidence = (v.confidence or "moderate").lower()

        if penetrance_score >= 0.8:
            # High penetrance (likely dominant): 1 copy sufficient
            zygosity_factor = min(kg.dosage, 1.0)
        elif penetrance_score >= 0.3:
            # Moderate penetrance: dosage-dependent
            zygosity_factor = min(kg.dosage * 0.7, 1.4)
        else:
            # Low/very_low penetrance
            if kg.dosage >= 2 and confidence in ("high", "moderate"):
                # Homozygous known risk variant with decent confidence.
                # 2.0× accounts for dose-dependent risk variants (e.g.
                # ABCG2 Q141K) that have strong effects when homozygous
                # despite low monogenic penetrance in the gene model.
                zygosity_factor = 2.0
            elif kg.dosage >= 2:
                zygosity_factor = 1.0
            else:
                zygosity_factor = 0.1  # heterozygous carrier → negligible

        contribution = v.contribution_score * zygosity_factor

        # ClinGen validity adjustment for Mendelian genes
        gene_w = gene_map.get(v.gene)
        if gene_w and gene_w.is_mendelian:
            clingen_mult = _clingen_weight(gene_w.clingen_validity, True)
            contribution *= clingen_mult

        results.append({
            "rsid": v.rsid or v.vcf_key,
            "gene": v.gene or "-",
            "variant": v.vcf_key,
            "gt": kg.gt,
            "dosage": kg.dosage,
            "risk_allele": v.effect_allele or v.alt,
            "contribution": round(min(contribution, 1.0), 3),
            "inferred_ref_ref": kg.inferred_ref_ref,
            "confidence": v.confidence,
            "note": v.note or "",
        })
    return results


def _score_dosage_risk(
    profile: DiseaseProfile,
    known_genotypes: list[KnownVariantGenotype],
) -> list[dict]:
    """Risk alleles with dosage effect."""
    results = []
    dosage_variants = [v for v in known_genotypes if v.variant.variant_class == "dosage_risk"]
    for kg in dosage_variants:
        v = kg.variant
        or_value = v.or_value or (math.exp(v.beta) if v.beta else None)
        if or_value and or_value > 0:
            per_allele_log_or = math.log(or_value)
        else:
            per_allele_log_or = 0.1
        contribution = kg.dosage * per_allele_log_or * v.contribution_score
        results.append({
            "rsid": v.rsid,
            "gene": v.gene,
            "variant": v.vcf_key,
            "gt": kg.gt,
            "risk_allele": v.effect_allele,
            "dosage": kg.dosage,
            "or_per_allele": or_value,
            "contribution": round(min(contribution, 1.0), 3),
            "inferred_ref_ref": kg.inferred_ref_ref,
            "note": v.note,
        })
    return results


def _score_gwas_prs(
    profile: DiseaseProfile,
    known_genotypes: list[KnownVariantGenotype],
    vcf_qc: dict,
) -> dict:
    """Weighted polygenic contribution from GWAS lead and PRS variants.

    Each variant contributes sqrt(|beta|) × dosage × contribution_score.
    Using sqrt(|beta|) amplifies small-effect variants (quantitative-trait
    GWAS with |beta| ~ 0.02) while preserving effect-direction ordering.
    The total is normalised by sqrt(n) for cross-disease comparability,
    with a VCF filter penalty applied when common variants are missing.
    """
    variants = [
        kg for kg in known_genotypes
        if kg.variant.variant_class in ("gwas_lead", "prs")
    ]
    if not variants:
        return {"score": 0.0, "percentile": None, "variant_count": 0, "variants": []}

    total = 0.0
    details = []
    for kg in variants:
        v = kg.variant
        beta = v.beta or (math.log(v.or_value) if v.or_value else 0.0)
        # Amplify small effect sizes via sqrt so quantitative-trait GWAS
        # (|beta| ~ 0.02) still produce meaningful per-locus contributions
        effective_beta = math.sqrt(abs(beta)) * (1 if beta >= 0 else -1)
        weight = v.contribution_score
        contribution = effective_beta * kg.dosage * weight
        total += contribution
        # VCF coordinates (may differ from template when REF/ALT are swapped)
        vcf_ref = kg.ref
        vcf_alt = kg.alt
        vcf_variant = f"{kg.chrom}:{kg.pos}:{vcf_ref}:{vcf_alt}"
        template_swap = (vcf_ref != v.ref or vcf_alt != v.alt) and not kg.inferred_ref_ref

        details.append({
            "rsid": v.rsid,
            "gene": v.gene,
            "variant": vcf_variant,  # use VCF coordinates for accurate display
            "template_variant": v.vcf_key if template_swap else None,
            "template_swap": template_swap,
            "gt": kg.gt,
            "effect_allele": v.effect_allele,
            "dosage": kg.dosage,
            "beta": beta,
            "beta_effective": round(effective_beta, 4),
            "weight": weight,
            "contribution": round(contribution, 4),
            "inferred_ref_ref": kg.inferred_ref_ref,
        })

    # Normalise by sqrt(variant count) for cross-disease comparability
    n = len(variants)
    total_normalised = total / math.sqrt(n)

    # Downweight if common variants are filtered in the input VCF
    if vcf_qc.get("likely_filtered"):
        total_normalised *= VCF_GWAS_DOWNWEIGHT_IF_FILTERED

    return {
        "score": round(total_normalised, 4),
        "score_raw": round(total, 4),
        "sqrt_n": round(math.sqrt(n), 1),
        "percentile": None,
        "variant_count": n,
        "variants": details,
    }


def _score_prs_high(
    profile: DiseaseProfile,
    known_genotypes: list[KnownVariantGenotype],
    vcf_qc: dict,
) -> dict:
    """High-confidence PRS variants with published weights.

    These are PRS variants from validated polygenic risk scores with
    independently replicated weights. We use sqrt(|beta|) amplification
    (same as gwas_prs) so that high-confidence small-effect variants still
    contribute meaningfully to the overall score, and return per-variant
    details for reporting.
    """
    variants = [
        kg for kg in known_genotypes
        if kg.variant.variant_class == "prs_high"
    ]
    if not variants:
        return {"score": 0.0, "variant_count": 0, "variants": []}

    total = 0.0
    details = []
    for kg in variants:
        v = kg.variant
        beta = v.beta or (math.log(v.or_value) if v.or_value else 0.0)
        # Use sqrt amplification to be consistent with gwas_prs; these are
        # well-replicated small-effect variants and should not be damped by
        # the raw beta linear scale.
        effective_beta = math.sqrt(abs(beta)) * (1 if beta >= 0 else -1)
        weight = v.contribution_score
        contribution = effective_beta * kg.dosage * weight
        total += contribution

        vcf_ref = kg.ref
        vcf_alt = kg.alt
        vcf_variant = f"{kg.chrom}:{kg.pos}:{vcf_ref}:{vcf_alt}"
        template_swap = (vcf_ref != v.ref or vcf_alt != v.alt) and not kg.inferred_ref_ref

        details.append({
            "rsid": v.rsid,
            "gene": v.gene,
            "variant": vcf_variant,
            "template_variant": v.vcf_key if template_swap else None,
            "template_swap": template_swap,
            "gt": kg.gt,
            "effect_allele": v.effect_allele,
            "dosage": kg.dosage,
            "beta": beta,
            "beta_effective": round(effective_beta, 4),
            "weight": weight,
            "contribution": round(contribution, 4),
            "inferred_ref_ref": kg.inferred_ref_ref,
            "note": v.note,
        })

    # Normalise by sqrt(n) for cross-disease comparability
    n = len(variants)
    total_normalised = total / math.sqrt(n)

    # Apply VCF filter penalty (lighter than gwas_prs because these are
    # high-confidence PRS variants whose weights have been independently
    # replicated).
    if vcf_qc.get("likely_filtered"):
        total_normalised *= VCF_GWAS_DOWNWEIGHT_IF_FILTERED

    return {
        "score": round(total_normalised, 4),
        "score_raw": round(total, 4),
        "sqrt_n": round(math.sqrt(n), 1),
        "variant_count": n,
        "variants": details,
    }


def _score_regulatory(
    profile: DiseaseProfile,
    tiered: dict[str, list[dict]],
) -> list[dict]:
    """Weak evidence: rare variants in regulatory regions of core genes."""
    results = []
    core_genes = profile.core_genes
    gene_map = {g.gene: g for g in profile.gene_set}
    for v in tiered.get("tier3_variants", []):
        gene = _variant_gene(v)
        if gene not in core_genes:
            continue
        if not _is_high_impact(v):
            continue
        af = _gnomad_af(v)
        if af is not None and af > 0.01:
            continue
        gw = gene_map.get(gene)
        contribution = (gw.contribution_score * 0.1) if gw else 0.05
        results.append({
            "gene": gene,
            "variant": _variant_key(v),
            "contribution": round(contribution, 3),
            "note": "Rare functional variant in regulatory region of core gene",
            "raw": v,
        })
    return results


def _collect_clinvar_enriched(
    tiered: dict[str, list[dict]],
) -> list[dict]:
    """Collect all variants that have ClinVar annotations for reporting.

    This layer is informational only and does not affect scoring.
    """
    results = []
    seen: set[str] = set()
    for tier_key in ("tier1_variants", "tier2_variants", "tier3_variants"):
        tier = 1 if tier_key == "tier1_variants" else (2 if tier_key == "tier2_variants" else 3)
        for v in tiered.get(tier_key, []):
            sig = v.get("clinvar_sig") or ""
            if not sig:
                continue
            key = v.get("variant_id") or f"{v.get('chrom')}:{v.get('pos')}:{v.get('ref')}:{v.get('alt')}"
            if key in seen:
                continue
            seen.add(key)
            results.append({
                "gene": _variant_gene(v),
                "variant": key,
                "tier": tier,
                "clinvar_sig": sig,
                "clinvar_category": v.get("clinvar_category", "other"),
                "clinvar_disease_match": v.get("clinvar_disease_match", False),
                "clinvar_diseases": v.get("clinvar_diseases", []),
                "note": "ClinVar annotated variant (informational)",
                "raw": v,
            })
    return results


def _layer_level(score: float, layer: str, profile: DiseaseProfile) -> str:
    """Return a per-layer qualitative level."""
    if score <= 0.0:
        return "none"
    # Mendelian layers are binary-presence driven
    if layer in ("mendelian_high", "mendelian_mod"):
        return "high" if score >= 0.8 else ("moderate" if score >= 0.3 else "low")
    # Known pathogenic: presence is strong evidence
    if layer == "known_pathogenic":
        return "high" if score >= 0.8 else ("moderate" if score >= 0.3 else "low")
    # Dosage and GWAS/PRS use continuous thresholds
    if layer == "dosage_risk":
        return "high" if score >= 0.5 else ("moderate" if score >= 0.15 else "low")
    if layer == "gwas_prs":
        return "high" if score >= 1.0 else ("moderate" if score >= 0.3 else "low")
    if layer == "regulatory":
        return "moderate" if score >= 0.1 else "low"
    return "low"


def _overall_level(score: Optional[float], profile: DiseaseProfile) -> str:
    if score is None:
        return "uncertain"
    if profile.mode == "mendelian":
        if score >= 0.8:
            return "high"
        if score >= 0.5:
            return "moderate"
        if score >= 0.2:
            return "low"
        return "very_low"
    else:
        # Complex mode: weights span 0–3.6, but in practice a single strong
        # layer (e.g. known_pathogenic homozygous + modifiers) reaches ~0.7.
        # Thresholds are set so that:
        #   high   = mendelian_high hit or multiple strong layers
        #   moderate = one strong known variant or clear polygenic signal
        #   low    = some genetic evidence but insufficient alone
        if score >= 1.0:
            return "high"
        if score >= 0.5:
            return "moderate"
        if score >= 0.2:
            return "low"
        return "very_low"


def score(
    profile: DiseaseProfile,
    tiered_variants: dict[str, list[dict]],
    known_genotypes: list[KnownVariantGenotype],
    vcf_qc: dict,
) -> ContributionResult:
    """Run the full contribution scoring model."""
    result = ContributionResult()

    result.mendelian_high = _score_mendelian_high(profile, tiered_variants)
    result.mendelian_mod = _score_mendelian_mod(profile, tiered_variants)
    result.known_pathogenic = _score_known_pathogenic(profile, known_genotypes)
    result.dosage_risk = _score_dosage_risk(profile, known_genotypes)
    result.gwas_prs = _score_gwas_prs(profile, known_genotypes, vcf_qc)
    result.regulatory = _score_regulatory(profile, tiered_variants)
    result.clinvar_enriched = _collect_clinvar_enriched(tiered_variants)

    # High-confidence PRS layer
    result.prs_high = _score_prs_high(profile, known_genotypes, vcf_qc)

    # Combine into overall score using model weights
    model = profile.contribution_model or {}
    high_score = sum(x["contribution"] for x in result.mendelian_high)
    mod_score = sum(x["contribution"] for x in result.mendelian_mod)
    known_score = sum(x["contribution"] for x in result.known_pathogenic)
    dosage_score = sum(x["contribution"] for x in result.dosage_risk)
    gwas_score = abs(result.gwas_prs.get("score", 0.0))
    reg_score = sum(x["contribution"] for x in result.regulatory)
    prs_high_score = abs(result.prs_high.get("score", 0.0))

    weights = {
        "mendelian_high": model.get("mendelian_high", {}).get("weight", 1.0),
        "mendelian_mod": model.get("mendelian_mod", {}).get("weight", 0.8),
        "known_pathogenic": model.get("known_pathogenic", {}).get("weight", 0.9),
        "prs_high": model.get("prs_high", {}).get("weight", 0.9),
        "dosage_risk": model.get("dosage_risk", {}).get("weight", 0.5),
        "gwas_prs": model.get("gwas_prs", {}).get("weight", 0.3),
        "regulatory": model.get("regulatory", {}).get("weight", 0.1),
    }

    overall = (
        min(high_score, 1.0) * weights["mendelian_high"]
        + min(mod_score, 1.0) * weights["mendelian_mod"]
        + min(known_score, 1.0) * weights["known_pathogenic"]
        + min(dosage_score, 1.0) * weights["dosage_risk"]
        + min(gwas_score, 1.0) * weights["gwas_prs"]
        + min(prs_high_score, 1.0) * weights["prs_high"]
        + min(reg_score, 1.0) * weights["regulatory"]
    )

    # If any high-penetrance variant found, floor the score
    if result.mendelian_high:
        overall = max(overall, 0.8)

    result.overall_score = round(overall, 3)
    result.overall_level = _overall_level(result.overall_score, profile)

    result.layer_levels = {
        "mendelian_high": _layer_level(high_score, "mendelian_high", profile),
        "mendelian_mod": _layer_level(mod_score, "mendelian_mod", profile),
        "known_pathogenic": _layer_level(known_score, "known_pathogenic", profile),
        "dosage_risk": _layer_level(dosage_score, "dosage_risk", profile),
        "gwas_prs": _layer_level(gwas_score, "gwas_prs", profile),
        "prs_high": _layer_level(prs_high_score, "prs_high", profile),
        "regulatory": _layer_level(reg_score, "regulatory", profile),
    }

    result.details = [
        {"layer": "mendelian_high", "count": len(result.mendelian_high), "score": round(high_score, 3)},
        {"layer": "mendelian_mod", "count": len(result.mendelian_mod), "score": round(mod_score, 3)},
        {"layer": "known_pathogenic", "count": len(result.known_pathogenic), "score": round(known_score, 3)},
        {"layer": "dosage_risk", "count": len(result.dosage_risk), "score": round(dosage_score, 3)},
        {"layer": "gwas_prs", "count": result.gwas_prs["variant_count"], "score": round(gwas_score, 3)},
        {"layer": "prs_high", "count": result.prs_high["variant_count"], "score": round(prs_high_score, 3)},
        {"layer": "regulatory", "count": len(result.regulatory), "score": round(reg_score, 3)},
        {"layer": "clinvar_enriched", "count": len(result.clinvar_enriched), "score": 0.0},
    ]

    return result
