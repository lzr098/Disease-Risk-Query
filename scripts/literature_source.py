"""Literature evidence source for disease risk queries.

Primary source: built-in key publications from constants.DISEASE_BUILTIN_REFS.
Optional extension: PubMed E-utilities queries can be added later.
"""

from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Optional

from constants import DISEASE_BUILTIN_REFS

logger = logging.getLogger(__name__)


def _normalize_disease_key(disease_name: str) -> str:
    return re.sub(r"[^a-z0-9 ]+", " ", disease_name.lower()).strip()


def get_key_literature(disease_name: str) -> list[dict]:
    """Return curated key literature entries for a disease."""
    key = _normalize_disease_key(disease_name)
    builtin = DISEASE_BUILTIN_REFS.get(key) or DISEASE_BUILTIN_REFS.get(disease_name.lower())
    if builtin:
        return builtin.get("key_literature", [])
    return []


def annotate_literature_support(
    variants: list[dict],
    disease_name: str,
) -> list[dict]:
    """Tag variants whose gene appears in key literature for the disease.

    Adds 'literature_support' (list of PMIDs/titles) to each variant.
    """
    lit = get_key_literature(disease_name)
    if not lit:
        for v in variants:
            v["literature_support"] = []
        return variants

    # Build gene -> evidence map
    gene_to_evidence: dict[str, list[dict]] = {}
    for entry in lit:
        for gene in entry.get("genes", []):
            gene_to_evidence.setdefault(gene.upper(), []).append(entry)

    for v in variants:
        gene = (v.get("gene") or v.get("GENE") or "").upper()
        v["literature_support"] = gene_to_evidence.get(gene, [])
    return variants


def score_literature_support(
    variants: list[dict],
    disease_name: str,
    core_genes: Optional[set[str]] = None,
) -> dict:
    """Summarize literature support across a variant list.

    Returns dict with:
      - variant_hits: variants with literature support
      - gene_hits: genes with literature support
      - core_genes_covered: overlap with core gene set
    """
    annotated = annotate_literature_support(variants, disease_name)
    hits = [v for v in annotated if v.get("literature_support")]
    genes = sorted({(v.get("gene") or v.get("GENE") or "").upper() for v in hits})
    core_overlap = sorted(set(genes) & (core_genes or set())) if core_genes else []
    return {
        "variant_hits": hits,
        "gene_hits": genes,
        "core_genes_covered": core_overlap,
        "total_entries": len(get_key_literature(disease_name)),
    }
