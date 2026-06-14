"""Literature evidence source for disease risk queries.

Primary source: built-in key publications from constants.DISEASE_BUILTIN_REFS.
Optional extension: PubMed E-utilities queries are used to supplement curated
literature when building a DiseaseProfile.
"""

from __future__ import annotations

import json
import logging
import re
import time
from pathlib import Path
from typing import Optional

from constants import DISEASE_BUILTIN_REFS

logger = logging.getLogger(__name__)

# NCBI E-utilities endpoints
PUBMED_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PUBMED_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


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


def _pubmed_query(
    disease_name: str,
    genes: list[str],
    retmax: int = 20,
) -> list[dict]:
    """Query PubMed E-utilities for disease+gene literature.

    Uses requests with trust_env=False to avoid proxy issues. Returns a list of
    literature entry dicts. On any error, returns an empty list so the pipeline
    can fall back to curated literature.
    """
    try:
        import requests
    except ImportError:
        logger.warning("requests not installed; cannot query PubMed")
        return []

    gene_terms = " OR ".join(f"{g}[Title/Abstract]" for g in genes[:10])
    term = f"({disease_name}[Title/Abstract]) AND ({gene_terms}) AND (variant OR mutation OR germline)"

    try:
        # ESearch
        search_resp = requests.get(
            PUBMED_ESEARCH,
            params={"db": "pubmed", "term": term, "retmax": retmax, "retmode": "json"},
            timeout=30,
        )
        search_resp.raise_for_status()
        search_data = search_resp.json()
        id_list = search_data.get("esearchresult", {}).get("idlist", [])
        if not id_list:
            return []

        # ESummary
        summary_resp = requests.get(
            PUBMED_ESUMMARY,
            params={"db": "pubmed", "id": ",".join(id_list), "retmode": "json"},
            timeout=30,
        )
        summary_resp.raise_for_status()
        summary_data = summary_resp.json()
        results = summary_data.get("result", {})

        entries = []
        for pmid in id_list:
            info = results.get(pmid, {})
            title = info.get("title", "")
            authors = [a.get("name", "") for a in info.get("authors", [])[:3]]
            source = info.get("source", "")
            pubdate = info.get("pubdate", "")
            year = pubdate.split()[0] if pubdate else ""
            entries.append({
                "pmid": str(pmid),
                "title": title,
                "authors": authors,
                "journal": source,
                "year": year,
                "genes": [],
                "note": f"PubMed search for {disease_name}",
                "source": "pubmed",
            })
        return entries
    except Exception as exc:
        logger.warning("PubMed query failed: %s", exc)
        return []


def build_literature_for_profile(
    disease_name: str,
    genes: list[str],
    cache_root: Optional[Path] = None,
    use_pubmed: bool = True,
    pubmax: int = 20,
) -> list[dict]:
    """Build the literature list for a DiseaseProfile.

    Combines curated built-in literature with a dynamic PubMed query. Results
    are cached under cache_root/literature/{disease}.json if cache_root is
    provided.

    Each entry has a 'source' field: 'curated' or 'pubmed'.
    """
    norm_key = _normalize_disease_key(disease_name)
    cached = []
    cache_path: Optional[Path] = None
    if cache_root:
        cache_path = Path(cache_root) / "literature" / f"{norm_key}.json"
        if cache_path.exists():
            try:
                cached = json.loads(cache_path.read_text(encoding="utf-8"))
                logger.info("Loaded %d cached literature entries for %s", len(cached), disease_name)
                return cached
            except Exception as exc:
                logger.warning("Failed to load literature cache: %s", exc)

    # Curated base
    curated = get_key_literature(disease_name)
    for entry in curated:
        entry.setdefault("source", "curated")
        entry.setdefault("genes", [])

    seen_pmids = {str(e.get("pmid")) for e in curated if e.get("pmid")}

    # Dynamic supplement
    dynamic: list[dict] = []
    if use_pubmed and genes:
        dynamic = _pubmed_query(disease_name, genes, retmax=pubmax)
        dynamic = [d for d in dynamic if d["pmid"] not in seen_pmids]
        # Assign overlapping genes (heuristic)
        gene_set = {g.upper() for g in genes}
        for d in dynamic:
            matched_genes = [g for g in gene_set if g.lower() in d["title"].lower()]
            d["genes"] = matched_genes[:3]

    combined = curated + dynamic

    if cache_path:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(combined, indent=2, ensure_ascii=False), encoding="utf-8")

    logger.info(
        "Built literature for %s: %d curated + %d dynamic = %d entries",
        disease_name, len(curated), len(dynamic), len(combined),
    )
    return combined
