"""Build or load a DiseaseProfile for a given disease query.

The builder supports two paths:
1. Built-in disease templates (constants.DISEASE_BUILTIN_REFS) are normalized
   into DiseaseProfile objects and cached on disk.
2. Diseases without a built-in template fall back to HPO + OMIM + ClinVar
   to build a minimal profile.
"""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
from typing import Any, Optional

from constants import (
    DISEASE_BUILTIN_REFS,
    DISEASE_REF_CACHE,
    resolve_builtin_disease_key,
    resolve_disease_mode,
)
from disease_profile import (
    DiseaseProfile,
    GeneDomainInfo,
    GeneWeight,
    RegulatoryRegion,
    VariantWeight,
    build_default_contribution_model,
)
from disease_reference import DiseaseReference
from gene_set_builder import GeneSetBuilder
from hpo_mapper import HPOMapper, resolve_disease_query
from literature_source import build_literature_for_profile

logger = logging.getLogger(__name__)


def _normalize_builtin_ref(key: str, ref: dict[str, Any]) -> dict[str, Any]:
    """Ensure a built-in disease reference has all DiseaseProfile fields.

    Backward-compatible: fills missing fields and converts legacy fields.
    """
    # Required top-level fields
    normalized: dict[str, Any] = {
        "canonical_name": key,
        "aliases": ref.get("aliases", []),
        "mode": ref.get("mode", resolve_disease_mode(key)),
        "hpo_terms": ref.get("hpo_terms", []),
        "gene_set": [],
        "known_pathogenic_variants": list(ref.get("known_pathogenic_variants", [])),
        "gwas_lead_snps": list(ref.get("gwas_lead_snps", [])),
        "prs_variants": list(ref.get("prs_variants", [])),
        "prs_variants_high": list(ref.get("prs_variants_high", [])),
        "regulatory_regions": [],
        "key_regions": {},
        "key_literature": list(ref.get("key_literature", [])),
        "contribution_model": ref.get("contribution_model", {}),
    }

    # gene_set: use new structured field or derive from legacy core/gwas genes
    gene_set = ref.get("gene_set", [])
    if gene_set:
        normalized["gene_set"] = list(gene_set)
    else:
        # Legacy fallback
        core = set(ref.get("core_genes", []))
        gwas = set(ref.get("gwas_loci_genes", []))
        normalized["gene_set"] = [
            {
                "gene": g,
                "tier": "mendelian_mod" if normalized["mode"] == "mendelian" else "mendelian_mod",
                "contribution_score": 1.0,
                "penetrance": ">0.95" if normalized["mode"] == "mendelian" else "moderate",
                "penetrance_score": 0.95 if normalized["mode"] == "mendelian" else 0.3,
                "evidence": "mixed",
                "note": "",
            }
            for g in sorted(core)
        ] + [
            {
                "gene": g,
                "tier": "gwas",
                "contribution_score": 0.15,
                "penetrance": "very_low",
                "penetrance_score": 0.05,
                "evidence": "gwas",
                "note": "",
            }
            for g in sorted(gwas - core)
        ]

    # regulatory_regions: map from legacy "regions" (genomic regions of interest)
    regions = ref.get("regions", [])
    normalized["regulatory_regions"] = [
        {
            "chrom": r.get("chrom", ""),
            "start": int(r.get("start", 0)),
            "end": int(r.get("end", 0)),
            "gene": r.get("gene") or None,
            "rtype": r.get("type", "regulatory"),
            "source": r.get("source", "builtin"),
            "note": r.get("note", ""),
        }
        for r in regions
        if r.get("chrom") and int(r.get("start", 0)) < int(r.get("end", 0))
    ]

    # key_regions: protein domain info for domain-dive
    key_regions = ref.get("key_regions", {})
    normalized["key_regions"] = {
        gene: {
            "note": info.get("note", ""),
            "regions": [
                {"name": r.get("name", ""), "residues": r.get("residues", ""), "note": r.get("note", "")}
                for r in info.get("regions", [])
            ],
            "critical_residues": [
                {"residue": r.get("residue", ""), "note": r.get("note", "")}
                for r in info.get("critical_residues", [])
            ],
        }
        for gene, info in key_regions.items()
    }

    # Ensure contribution_model exists
    if not normalized["contribution_model"]:
        normalized["contribution_model"] = build_default_contribution_model(normalized["mode"])

    return normalized


def _dict_to_disease_profile(data: dict[str, Any]) -> DiseaseProfile:
    """Convert a normalized dict into a DiseaseProfile dataclass."""
    return DiseaseProfile.from_dict(data)


def _enrich_with_hpo_omim(
    profile: DiseaseProfile,
    disease_name: str,
    hpo_mapper: Optional[HPOMapper] = None,
    gene_builder: Optional[GeneSetBuilder] = None,
    max_extra_genes: int = 50,
) -> DiseaseProfile:
    """Add HPO/OMIM-derived genes if the built-in template is sparse."""
    hpo_mapper = hpo_mapper or HPOMapper()
    gene_builder = gene_builder or GeneSetBuilder()

    hpo_result = resolve_disease_query(disease_name)
    hpo_id = hpo_result.get("hpo_id")
    hpo_genes = []
    if hpo_id:
        hpo_genes = sorted(hpo_mapper.genes_for_hpo(hpo_id))

    canonical = resolve_builtin_disease_key(disease_name) or disease_name
    omim_genes = gene_builder.query_omim([canonical, disease_name], max_genes=max_extra_genes * 2)

    existing = profile.all_genes
    new_genes = (set(hpo_genes) | omim_genes) - existing
    if not new_genes:
        return profile

    logger.info("Enriching profile with %d HPO/OMIM genes", len(new_genes))
    for gene in sorted(new_genes):
        profile.gene_set.append(
            GeneWeight(
                gene=gene,
                tier="gwas",
                contribution_score=0.15,
                penetrance="very_low",
                penetrance_score=0.05,
                evidence="hpo_omim",
                note="Added from HPO/OMIM lookup",
            )
        )
    return profile


def build_or_load_profile(
    disease_name: str,
    hpo_id: Optional[str] = None,
    refresh: bool = False,
    cache_root: Path = DISEASE_REF_CACHE,
    enrich: Optional[bool] = None,
) -> DiseaseProfile:
    """Return a DiseaseProfile for the disease query.

    Uses built-in template if available; otherwise builds a minimal profile
    from HPO + OMIM + ClinVar.
    """
    canonical = resolve_builtin_disease_key(disease_name)
    cache_key = canonical or disease_name.lower().replace(" ", "_").replace("/", "_")
    cache_dir = cache_root / cache_key
    profile_path = cache_dir / "profile.json"

    if profile_path.exists() and not refresh:
        # Auto-invalidate cache when the built-in template has changed
        if canonical and canonical in DISEASE_BUILTIN_REFS:
            current_hash = hashlib.md5(
                json.dumps(DISEASE_BUILTIN_REFS[canonical], sort_keys=True).encode()
            ).hexdigest()
            hash_path = cache_dir / "profile.template_hash"
            cached_hash = hash_path.read_text().strip() if hash_path.exists() else None
            if cached_hash != current_hash:
                logger.info(
                    "Built-in template changed for '%s'; rebuilding profile "
                    "(cached hash %s → current %s)", canonical, cached_hash, current_hash
                )
                hash_path.write_text(current_hash)
                # Fall through to rebuild
            else:
                logger.info("Loading cached DiseaseProfile from %s", profile_path)
                return DiseaseProfile.from_json(profile_path)
        else:
            logger.info("Loading cached DiseaseProfile from %s", profile_path)
            return DiseaseProfile.from_json(profile_path)

    profile: DiseaseProfile
    if canonical and canonical in DISEASE_BUILTIN_REFS:
        logger.info("Building DiseaseProfile from built-in template for '%s'", disease_name)
        normalized = _normalize_builtin_ref(canonical, DISEASE_BUILTIN_REFS[canonical])
        profile = _dict_to_disease_profile(normalized)
        if hpo_id and not profile.hpo_terms:
            profile.hpo_terms.append({"id": hpo_id, "name": disease_name})
    else:
        logger.info("No built-in template for '%s'; building minimal profile", disease_name)
        mode = resolve_disease_mode(disease_name)
        profile = DiseaseProfile(
            canonical_name=disease_name,
            mode=mode,
            contribution_model=build_default_contribution_model(mode),
        )
        if hpo_id:
            profile.hpo_terms.append({"id": hpo_id, "name": disease_name})

    # Default enrichment behavior: enrich only when there is no built-in template,
    # so built-in templates stay clean and fast.
    canonical = resolve_builtin_disease_key(disease_name)
    should_enrich = enrich if enrich is not None else (canonical is None)
    if should_enrich:
        profile = _enrich_with_hpo_omim(profile, disease_name)

    # Build literature evidence for the profile (curated + dynamic PubMed).
    # This runs for both built-in templates and custom-built profiles.
    try:
        profile.key_literature = build_literature_for_profile(
            disease_name=disease_name,
            genes=sorted(profile.all_genes),
            cache_root=cache_root,
            use_pubmed=True,
            pubmax=20,
        )
    except Exception as exc:
        logger.warning("Failed to build literature for profile: %s", exc)

    # Backward-compatible disease reference cache is built lazily by callers that
    # still use get_disease_reference(); we avoid forcing it here to keep profile
    # construction fast.

    cache_dir.mkdir(parents=True, exist_ok=True)
    profile.to_json(profile_path)
    logger.info("Saved DiseaseProfile to %s", profile_path)

    # Save template hash for cache invalidation on next load
    if canonical and canonical in DISEASE_BUILTIN_REFS:
        current_hash = hashlib.md5(
            json.dumps(DISEASE_BUILTIN_REFS[canonical], sort_keys=True).encode()
        ).hexdigest()
        (cache_dir / "profile.template_hash").write_text(current_hash)

    return profile


def list_builtin_profiles() -> list[str]:
    """Return list of diseases with built-in templates."""
    return sorted(DISEASE_BUILTIN_REFS.keys())
