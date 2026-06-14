"""Post-GPA tier filtering and denoising."""

from __future__ import annotations

import logging
from collections import defaultdict
from typing import Optional

from constants import (
    AD_CORE_GENES,
    TIER3_CORE_GENE_MAX_GNOMAD_AF,
    TIER3_KEEP_IMPACTS,
    TIER3_MAX_GNOMAD_AF,
    TIER3_MAX_VARIANTS_PER_GENE,
)

# Stricter thresholds for GWAS/association loci (common variants are expected)
TIER3_GWAS_LOCUS_MAX_GNOMAD_AF = 0.005
TIER3_GWAS_LOCUS_MAX_VARIANTS = 5

logger = logging.getLogger(__name__)


def _gnomad_af(variant: dict) -> Optional[float]:
    """Extract gnomAD AF from a GPA variant dict."""
    for key in ("gnomad_af", "gnomAD_AF", "AF", "af"):
        val = variant.get(key)
        if val is not None and val != "":
            try:
                return float(val)
            except (TypeError, ValueError):
                continue
    return None


def _impact(variant: dict) -> str:
    """Return IMPACT/impact field uppercase."""
    for key in ("IMPACT", "impact", "primary_impact"):
        val = variant.get(key)
        if val:
            return str(val).upper()
    return ""


def _gene(variant: dict) -> str:
    for key in ("GENE", "gene"):
        val = variant.get(key)
        if val:
            return str(val)
    return ""


def _is_likely_benign_clinvar(variant: dict) -> bool:
    """Return True if ClinVar annotation is explicitly benign without pathogenic evidence."""
    import re
    sig = str(variant.get("clinvar") or variant.get("CLIN_SIG") or "").lower()
    if not sig:
        return False
    # If it contains pathogenic or conflicting, do not treat as benign
    if any(x in sig for x in ("pathogenic", "conflicting", "drug_response", "risk_factor")):
        return False
    return "benign" in sig or "likely_benign" in sig


def denoise_tier3(
    tier3_variants: list[dict],
    core_genes: Optional[set[str]] = None,
    gwas_loci_genes: Optional[set[str]] = None,
) -> list[dict]:
    """Denoise Tier 3 by AF, impact, ClinVar, and per-gene caps.

    Rules:
      - Keep only HIGH/MODERATE impact variants.
      - Drop variants explicitly annotated as benign/likely_benign unless the
        gene is a Mendelian core gene (where clinical interpretation is more
        nuanced).
      - For Mendelian core genes, allow AF up to TIER3_CORE_GENE_MAX_GNOMAD_AF.
      - For GWAS/association loci, apply TIER3_GWAS_LOCUS_MAX_GNOMAD_AF and a
        lower per-gene cap because these regions are often highly polymorphic
        (e.g. HLA) and we are mainly interested in rare damaging variants.
      - For other genes, require AF < TIER3_MAX_GNOMAD_AF.
      - Cap each gene to the configured maximum (rarest first).
    """
    core_genes = core_genes or AD_CORE_GENES
    gwas_loci_genes = gwas_loci_genes or set()
    kept = []
    for v in tier3_variants:
        impact = _impact(v)
        if impact not in TIER3_KEEP_IMPACTS:
            continue
        gene = _gene(v)
        # Drop clearly benign variants outside core genes
        if gene not in core_genes and _is_likely_benign_clinvar(v):
            continue
        af = _gnomad_af(v)
        if gene in core_genes:
            max_af = TIER3_CORE_GENE_MAX_GNOMAD_AF
        elif gene in gwas_loci_genes:
            max_af = TIER3_GWAS_LOCUS_MAX_GNOMAD_AF
        else:
            max_af = TIER3_MAX_GNOMAD_AF
        if af is None:
            # Keep if no AF available and impact is high/moderate
            kept.append(v)
        elif af <= max_af:
            kept.append(v)

    # Per-gene cap: keep rarest variants first
    by_gene: defaultdict[str, list[tuple[Optional[float], dict]]] = defaultdict(list)
    for v in kept:
        by_gene[_gene(v)].append((_gnomad_af(v), v))

    final = []
    for gene, items in by_gene.items():
        # Sort by AF ascending (None treated as 0, i.e. rarest)
        items.sort(key=lambda x: x[0] if x[0] is not None else -1.0)
        if gene in gwas_loci_genes:
            cap = TIER3_GWAS_LOCUS_MAX_VARIANTS
        else:
            cap = TIER3_MAX_VARIANTS_PER_GENE
        if len(items) > cap:
            logger.info(
                "Tier 3 denoising: capped %s from %d to %d rarest variants",
                gene, len(items), cap,
            )
        for _, v in items[:cap]:
            final.append(v)

    logger.info(
        "Tier 3 denoising: %d -> %d variants (%d genes)",
        len(tier3_variants), len(final), len(by_gene),
    )
    return final


def _variant_key(v: dict) -> str:
    chrom = v.get("chrom") or v.get("CHROM") or ""
    pos = v.get("pos") or v.get("POS") or ""
    ref = v.get("ref") or v.get("REF") or ""
    alt = v.get("alt") or v.get("ALT") or ""
    return f"{chrom}:{pos}:{ref}:{alt}"


def remove_duplicate_variants(variants: list[dict]) -> list[dict]:
    """Remove duplicate variants by chrom:pos:ref:alt, keeping first occurrence."""
    seen: set[str] = set()
    out = []
    for v in variants:
        key = _variant_key(v)
        if key in seen:
            continue
        seen.add(key)
        out.append(v)
    return out
