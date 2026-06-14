"""Build disease-associated gene set from OMIM + HPO + optional HGNC validation."""

from __future__ import annotations

import json
import logging
import re
import sqlite3
import time
from pathlib import Path
from typing import Iterable, Optional

from constants import AD_CORE_GENES, HGNC_LOOKUP, OMIM_DB, OMIM_TITLE_GENE_CACHE

logger = logging.getLogger(__name__)


def _safe_split_symbols(text: str) -> Iterable[str]:
    """Yield cleaned gene-like tokens from an OMIM field."""
    if not text:
        return
    for sym in re.split(r"[,;|]", str(text)):
        sym = sym.strip()
        if sym:
            yield sym


class GeneSetBuilder:
    """Assemble and standardize a disease gene set."""

    def __init__(
        self,
        omim_db: Path = OMIM_DB,
        hgnc_lookup: Path = HGNC_LOOKUP,
    ) -> None:
        self.omim_db = omim_db
        self.hgnc_lookup = hgnc_lookup
        self._hgnc: dict[str, str] = {}
        self._omim_cache: Optional[list[dict]] = None
        self._omim_cache_loaded_at: float = 0.0

    def _load_hgnc(self) -> None:
        if self._hgnc:
            return
        if not self.hgnc_lookup.exists():
            logger.warning("HGNC lookup not found: %s", self.hgnc_lookup)
            return
        try:
            with open(self.hgnc_lookup, "r", encoding="utf-8") as f:
                data = json.load(f)
            # The local HGNC file is {"metadata": ..., "lookup": {"SYM": {...}, ...}}
            if isinstance(data, dict) and "lookup" in data:
                data = data["lookup"]
            # Expect {symbol -> approved_symbol} or nested dict
            if isinstance(data, dict):
                for k, v in data.items():
                    if isinstance(v, str):
                        self._hgnc[k.upper()] = v
                    elif isinstance(v, dict):
                        approved = v.get("approved_symbol") or v.get("symbol")
                        if approved:
                            self._hgnc[k.upper()] = approved
        except Exception as exc:
            logger.warning("Failed to load HGNC lookup: %s", exc)

    def standardize_symbol(self, symbol: str) -> str:
        """Return HGNC-approved symbol if known, else original."""
        self._load_hgnc()
        if not self._hgnc:
            return symbol
        return self._hgnc.get(symbol.upper(), symbol)

    def _load_omim_cache(self) -> list[dict]:
        """Load a compact OMIM title/symbol cache, rebuilding if stale or missing.

        The full OMIM table contains many large TEXT columns. Scanning it for every
        disease query is very slow. This cache keeps only the small fields needed
        for keyword matching (title, gene symbols, geneMap, phenotypeMap) and is
        rebuilt automatically when the source OMIM database changes.
        """
        if self._omim_cache is not None:
            db_mtime = self.omim_db.stat().st_mtime if self.omim_db.exists() else 0
            if self._omim_cache_loaded_at >= db_mtime:
                return self._omim_cache

        if (
            OMIM_TITLE_GENE_CACHE.exists()
            and self.omim_db.exists()
            and OMIM_TITLE_GENE_CACHE.stat().st_mtime >= self.omim_db.stat().st_mtime
        ):
            try:
                with open(OMIM_TITLE_GENE_CACHE, "r", encoding="utf-8") as f:
                    self._omim_cache = json.load(f)
                self._omim_cache_loaded_at = OMIM_TITLE_GENE_CACHE.stat().st_mtime
                logger.info("Loaded OMIM title/gene cache (%d entries)", len(self._omim_cache))
                return self._omim_cache
            except Exception as exc:
                logger.warning("Failed to load OMIM cache, rebuilding: %s", exc)

        logger.info("Building OMIM title/gene cache from %s", self.omim_db)
        cache: list[dict] = []
        if not self.omim_db.exists():
            logger.warning("OMIM DB not found: %s", self.omim_db)
            self._omim_cache = cache
            self._omim_cache_loaded_at = time.time()
            return cache

        try:
            conn = sqlite3.connect(str(self.omim_db))
            cur = conn.cursor()
            cur.execute(
                "SELECT title, hgnc_gene_symbol, hgnc_approved_gene_symbol, symbols, "
                "       geneMap, phenotypeMap "
                "FROM omim"
            )
            for title, hgnc, approved, symbols, gene_map, phenotype_map in cur:
                cache.append({
                    "title": title or "",
                    "hgnc_gene_symbol": hgnc or "",
                    "hgnc_approved_gene_symbol": approved or "",
                    "symbols": symbols or "",
                    "geneMap": gene_map or "",
                    "phenotypeMap": phenotype_map or "",
                })
            conn.close()

            OMIM_TITLE_GENE_CACHE.parent.mkdir(parents=True, exist_ok=True)
            with open(OMIM_TITLE_GENE_CACHE, "w", encoding="utf-8") as f:
                json.dump(cache, f)
            self._omim_cache = cache
            self._omim_cache_loaded_at = OMIM_TITLE_GENE_CACHE.stat().st_mtime
            logger.info("Built OMIM title/gene cache (%d entries)", len(cache))
        except Exception as exc:
            logger.warning("Failed to build OMIM cache: %s", exc)
            self._omim_cache = cache
            self._omim_cache_loaded_at = time.time()

        return cache

    def query_omim(
        self,
        keywords: Iterable[str],
        max_genes: int = 200,
    ) -> set[str]:
        """Search OMIM title/gene fields for disease keywords and return gene symbols."""
        genes: set[str] = set()
        cache = self._load_omim_cache()
        if not cache:
            logger.warning("OMIM cache empty, cannot search OMIM")
            return genes

        for kw in keywords:
            kw_clean = kw.strip().lower()
            if not kw_clean:
                continue
            for entry in cache:
                # Fast substring search on the small text fields
                if (
                    kw_clean in entry["title"].lower()
                    or kw_clean in entry["hgnc_gene_symbol"].lower()
                    or kw_clean in entry["hgnc_approved_gene_symbol"].lower()
                    or kw_clean in entry["symbols"].lower()
                    or kw_clean in entry["geneMap"].lower()
                    or kw_clean in entry["phenotypeMap"].lower()
                ):
                    for field in (
                        entry["hgnc_gene_symbol"],
                        entry["hgnc_approved_gene_symbol"],
                        entry["symbols"],
                    ):
                        for sym in _safe_split_symbols(field):
                            genes.add(sym)
                    if len(genes) >= max_genes * 4:
                        break

        return genes

    def build_gene_set(
        self,
        disease_name: str,
        hpo_genes: Optional[list[str]] = None,
        omim_keywords: Optional[list[str]] = None,
        extra_genes: Optional[list[str]] = None,
        max_genes: int = 200,
        core_genes: Optional[set[str]] = None,
    ) -> dict:
        """Build merged, standardized disease gene set.

        Args:
            disease_name: Primary disease name (also used as OMIM keyword).
            hpo_genes: Genes associated with the target HPO term.
            omim_keywords: Additional keywords for OMIM search.
            extra_genes: Literature-derived or manually provided genes.
            max_genes: Hard cap on returned gene count (safety limit).
            core_genes: Optional set of disease-specific core genes to prioritize.

        Returns:
            Dict with sources, gene lists, and final merged set.
        """
        hpo_genes = hpo_genes or []
        omim_keywords = list(omim_keywords or []) + [disease_name]
        extra_genes = extra_genes or []
        core_genes = core_genes or set()

        # OMIM
        omim_raw = self.query_omim(omim_keywords, max_genes=max_genes * 4)
        omim_genes = sorted(omim_raw)

        # Standardize all sources
        standardized_hpo = sorted({self.standardize_symbol(g) for g in hpo_genes})
        standardized_omim = sorted({self.standardize_symbol(g) for g in omim_genes})
        standardized_extra = sorted({self.standardize_symbol(g) for g in extra_genes})

        hpo_set = set(standardized_hpo)
        omim_set = set(standardized_omim)
        extra_set = set(standardized_extra)
        merged = hpo_set | omim_set | extra_set

        # Remove common false positives / non-gene strings
        non_gene_blacklist = {
            "ALZHEIMER", "DISEASE", "TYPE", "FORM", "FAMILIAL", "SPORADIC",
            "AMYLOIDOSIS", "CEREBRAL", "DEGENERATIVE", "PROTEIN", "GENE",
            "CHROMOSOME", "LINKED", "ASSOCIATED", "RELATED", "SUSCEPTIBILITY",
        }
        merged = {
            g for g in merged
            if re.match(r"^[A-Za-z0-9\-_\.]+$", g)
            and len(g) >= 2
            and not g.isdigit()
            and not any(word in g.upper() for word in non_gene_blacklist)
            and (not self._hgnc or g.upper() in self._hgnc)
        }

        # Rank genes by relevance instead of arbitrary alphabetical truncation.
        def _gene_priority(gene: str) -> tuple:
            score = 0
            if gene in core_genes:
                score += 100
            if gene in extra_set:
                score += 40
            if gene in omim_set:
                score += 20
            if gene in hpo_set:
                score += 10
            return (-score, gene)

        merged_sorted = sorted(merged, key=_gene_priority)
        truncated = False
        if len(merged_sorted) > max_genes:
            logger.warning(
                "Gene set truncated from %d to %d by --max-genes; highest-priority genes retained.",
                len(merged_sorted), max_genes,
            )
            merged_sorted = merged_sorted[:max_genes]
            truncated = True

        return {
            "disease_name": disease_name,
            "hpo_genes": standardized_hpo,
            "omim_genes": standardized_omim,
            "extra_genes": standardized_extra,
            "merged_genes": merged_sorted,
            "total": len(merged_sorted),
            "truncated": truncated,
            "sources": {
                "hpo": len(standardized_hpo),
                "omim": len(standardized_omim),
                "extra": len(standardized_extra),
            },
        }


def build_disease_gene_set(
    disease_name: str,
    hpo_genes: Optional[list[str]] = None,
    omim_keywords: Optional[list[str]] = None,
    extra_genes: Optional[list[str]] = None,
    max_genes: int = 200,
    core_genes: Optional[set[str]] = None,
) -> dict:
    """Convenience wrapper."""
    builder = GeneSetBuilder()
    return builder.build_gene_set(
        disease_name=disease_name,
        hpo_genes=hpo_genes,
        omim_keywords=omim_keywords,
        extra_genes=extra_genes,
        max_genes=max_genes,
        core_genes=core_genes,
    )
