"""Disease / phenotype query to HPO ID and gene set mapping.

This module does NOT rely on LLM; it performs deterministic mapping using:
1. Explicit HPO IDs passed by caller.
2. A curated disease-name -> HPO ID map for common conditions.
3. Keyword matching against HPO phenotype names from genes_to_phenotype.txt.

The calling agent is responsible for any natural-language disambiguation
before invoking this skill.
"""

from __future__ import annotations

import csv
import json
import logging
import re
from collections import defaultdict
from pathlib import Path
from typing import Optional

from constants import HPO_GENES_TO_PHENOTYPE, HPO_DISEASE_CACHE, resolve_builtin_disease_key

logger = logging.getLogger(__name__)


# Curated common disease -> HPO ID map (expand as needed)
DISEASE_TO_HPO: dict[str, str] = {
    "alzheimer disease": "HP:0000726",
    "alzheimer": "HP:0000726",
    "parkinson disease": "HP:0001300",
    "parkinson": "HP:0001300",
    "epilepsy": "HP:0001250",
    "seizure": "HP:0001250",
    "cardiomyopathy": "HP:0001638",
    "dilated cardiomyopathy": "HP:0001644",
    "hypertrophic cardiomyopathy": "HP:0001639",
    "arrhythmia": "HP:0011675",
    "long qt syndrome": "HP:0001657",
    "polycystic kidney disease": "HP:0000107",
    "nephrotic syndrome": "HP:0000100",
    "intellectual disability": "HP:0001249",
    "developmental delay": "HP:0001263",
    "autism": "HP:0000717",
    "schizophrenia": "HP:0100753",
    "bipolar disorder": "HP:0100753",
    "hereditary breast ovarian cancer": "HP:0003002",
    "breast cancer": "HP:0003002",
    "lynch syndrome": "HP:0003003",
    "colorectal cancer": "HP:0003003",
    "familial hypercholesterolemia": "HP:0003124",
    "marfan syndrome": "HP:0001654",
    "aortic aneurysm": "HP:0004942",
    "sickle cell disease": "HP:0001903",
    "thalassemia": "HP:0001903",
    "hemophilia": "HP:0001873",
    "type 1 diabetes": "HP:0000819",
    "type 2 diabetes": "HP:0000819",
    "maturity onset diabetes of the young": "HP:0000819",
    "glaucoma": "HP:0000501",
    "retinitis pigmentosa": "HP:0000510",
    "hearing loss": "HP:0000365",
    "deafness": "HP:0000365",
    "muscular dystrophy": "HP:0003560",
    "duchenne muscular dystrophy": "HP:0003560",
    "spinal muscular atrophy": "HP:0003202",
    "amyotrophic lateral sclerosis": "HP:0007354",
    "als": "HP:0007354",
    "huntington disease": "HP:0002185",
    "ataxia": "HP:0001251",
    "neuropathy": "HP:0009830",
    "migraine": "HP:0002076",
    "coagulation disorder": "HP:0001928",
    "bleeding disorder": "HP:0001928",
    "anemia": "HP:0001903",
    "leukemia": "HP:0001909",
    "lymphoma": "HP:0002665",
    "wilson disease": "HP:0001392",
    "hemochromatosis": "HP:0003272",
    "alpha-1 antitrypsin deficiency": "HP:0002032",
    "phenylketonuria": "HP:0001250",
    "hyperuricemia": "HP:0002149",
    "gout": "HP:0001997",
    "高尿酸血症": "HP:0002149",
    "痛风": "HP:0001997",
    "myocardial infarction": "HP:0001658",
    "mi": "HP:0001658",
    "coronary artery disease": "HP:0001658",
    "cad": "HP:0001658",
    "coronary heart disease": "HP:0001658",
    "chd": "HP:0001658",
    # Body mass / obesity continuous-trait mapping
    "body mass index": "HP:0045081",
    "bmi": "HP:0045081",
    "obesity": "HP:0001513",
    "obese": "HP:0001513",
    "肥胖": "HP:0001513",
    "肥胖倾向": "HP:0001513",
}


def _normalize(text: str) -> str:
    return re.sub(r"[^a-z0-9 ]+", " ", text.lower()).strip()


class HPOMapper:
    """Maps disease queries to HPO IDs and associated genes."""

    def __init__(self, hpo_file: Path = HPO_GENES_TO_PHENOTYPE) -> None:
        self.hpo_file = hpo_file
        self._hpo_to_genes: dict[str, set[str]] = defaultdict(set)
        self._hpo_to_name: dict[str, str] = {}
        self._loaded = False
        self._gene_cache: dict[str, list[str]] = {}

    def _load(self) -> None:
        if self._loaded:
            return
        if not self.hpo_file.exists():
            logger.warning("HPO file not found: %s", self.hpo_file)
            self._loaded = True
            return

        logger.info("Loading HPO gene-phenotype mapping from %s", self.hpo_file)
        with open(self.hpo_file, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                hpo_id = row.get("hpo_id", "").strip()
                hpo_name = row.get("hpo_name", "").strip()
                gene = row.get("gene_symbol", "").strip()
                if not hpo_id:
                    continue
                if gene:
                    self._hpo_to_genes[hpo_id].add(gene)
                if hpo_name and hpo_id not in self._hpo_to_name:
                    self._hpo_to_name[hpo_id] = hpo_name
        self._loaded = True

    def resolve(
        self,
        query: str,
        hpo_id: Optional[str] = None,
    ) -> dict:
        """Resolve a disease query to HPO ID, name, and gene set.

        Args:
            query: Natural-language disease name (e.g. "Alzheimer disease").
            hpo_id: Optional explicit HPO ID to use instead of keyword search.

        Returns:
            Dict with hpo_id, hpo_name, query, matched_by, genes.
        """
        self._load()

        if hpo_id:
            hpo_id = hpo_id.strip().upper()
            hpo_name = self._hpo_to_name.get(hpo_id, "Unknown HPO term")
            genes = sorted(self._hpo_to_genes.get(hpo_id, set()))
            return {
                "hpo_id": hpo_id,
                "hpo_name": hpo_name,
                "query": query,
                "matched_by": "explicit",
                "genes": genes,
            }

        # 0. Built-in disease alias resolution (handles Chinese names)
        canonical = resolve_builtin_disease_key(query)
        if canonical and canonical in DISEASE_TO_HPO:
            hpo_id = DISEASE_TO_HPO[canonical]
            return {
                "hpo_id": hpo_id,
                "hpo_name": self._hpo_to_name.get(hpo_id, "Unknown HPO term"),
                "query": query,
                "matched_by": "curated_builtin",
                "genes": sorted(self._hpo_to_genes.get(hpo_id, set())),
            }

        norm = _normalize(query)

        # 1. Curated direct map
        if norm in DISEASE_TO_HPO:
            hpo_id = DISEASE_TO_HPO[norm]
            return {
                "hpo_id": hpo_id,
                "hpo_name": self._hpo_to_name.get(hpo_id, "Unknown HPO term"),
                "query": query,
                "matched_by": "curated",
                "genes": sorted(self._hpo_to_genes.get(hpo_id, set())),
            }

        # 2. Curated keyword partial match (skip empty norm to avoid false matches)
        if norm:
            for disease_name, candidate_hpo in DISEASE_TO_HPO.items():
                if disease_name in norm or norm in disease_name:
                    return {
                        "hpo_id": candidate_hpo,
                        "hpo_name": self._hpo_to_name.get(candidate_hpo, "Unknown HPO term"),
                        "query": query,
                        "matched_by": "curated_partial",
                        "genes": sorted(self._hpo_to_genes.get(candidate_hpo, set())),
                    }

        # 3. Search HPO phenotype names for keyword overlap
        query_tokens = set(norm.split())
        best_hpo: Optional[str] = None
        best_score = 0
        for hid, name in self._hpo_to_name.items():
            name_tokens = set(_normalize(name).split())
            if not name_tokens:
                continue
            overlap = len(query_tokens & name_tokens)
            score = overlap / max(len(query_tokens), len(name_tokens))
            if score > best_score and score >= 0.5:
                best_score = score
                best_hpo = hid

        if best_hpo:
            return {
                "hpo_id": best_hpo,
                "hpo_name": self._hpo_to_name.get(best_hpo, ""),
                "query": query,
                "matched_by": "hpo_keyword",
                "match_score": round(best_score, 3),
                "genes": sorted(self._hpo_to_genes.get(best_hpo, set())),
            }

        return {
            "hpo_id": None,
            "hpo_name": None,
            "query": query,
            "matched_by": "none",
            "genes": [],
        }

    def genes_for_hpo(self, hpo_id: str) -> list[str]:
        self._load()
        hid = hpo_id.upper()
        if hid in self._gene_cache:
            return self._gene_cache[hid]
        result = sorted(self._hpo_to_genes.get(hid, set()))
        self._gene_cache[hid] = result
        return result

    def save_curated_cache(self, path: Path = HPO_DISEASE_CACHE) -> None:
        """Persist curated disease->HPO map for offline use."""
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w", encoding="utf-8") as f:
            json.dump(DISEASE_TO_HPO, f, indent=2, ensure_ascii=False)


_MAPPER_INSTANCE: Optional[HPOMapper] = None


def resolve_disease_query(query: str, hpo_id: Optional[str] = None) -> dict:
    """Convenience wrapper — reuses a module-level singleton."""
    global _MAPPER_INSTANCE
    if _MAPPER_INSTANCE is None:
        _MAPPER_INSTANCE = HPOMapper()
    return _MAPPER_INSTANCE.resolve(query, hpo_id)
