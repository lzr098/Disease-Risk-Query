"""Disease reference cache manager.

For each disease query, we build (or reuse) a persistent reference set that
aggregates:
  - core genes (familial + GWAS loci)
  - known pathogenic / likely pathogenic variants from ClinVar with phenotype match
  - GWAS lead SNPs (from built-in seeds or optional external query)
  - key literature evidence

These references are stored under ~/.workbuddy/data/drq_disease_refs/{disease_key}/
and reused across runs. New disease queries first check the cache; if absent,
the cache is built on demand.
"""

from __future__ import annotations

import gzip
import json
import logging
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from constants import (
    CLINVAR_VCF,
    DISEASE_BUILTIN_REFS,
    DISEASE_REF_CACHE,
    GENE_COORDS_CACHE,
    resolve_builtin_disease_key,
)
from clinvar_phenotype_matcher import disease_keywords

logger = logging.getLogger(__name__)


def _normalize_disease_key(disease_name: str) -> str:
    """Convert a disease query to a filesystem-safe cache key.

    Preserves non-ASCII characters (e.g., Chinese disease names) while
    removing filesystem-unsafe symbols.
    """
    # Replace filesystem-unsafe characters with underscore
    norm = re.sub(r'[\\/:*?"<>|]+', "_", disease_name.lower())
    # Normalize whitespace to single underscore
    norm = re.sub(r"\s+", "_", norm.strip())
    # Remove leading/trailing underscores
    norm = norm.strip("_")
    return norm


def _gene_to_bed_records(gene: str, coords_cache: Path = GENE_COORDS_CACHE) -> list[str]:
    """Fetch BED lines for a gene symbol from the GENCODE gene coords cache."""
    if not coords_cache.exists():
        return []
    try:
        data = json.loads(coords_cache.read_text(encoding="utf-8"))
    except Exception:
        return []
    records = []
    for chrom, start, end, strand in data.get(gene, []):
        records.append(f"{chrom}\t{start}\t{end}\t{gene}\t{strand}")
    return records


def _build_genes_bed(genes: list[str], output: Path) -> int:
    """Write a BED of gene loci for the reference set."""
    output.parent.mkdir(parents=True, exist_ok=True)
    seen: set[str] = set()
    count = 0
    with open(output, "w", encoding="utf-8") as f:
        for gene in sorted(set(genes)):
            for line in _gene_to_bed_records(gene):
                key = "\t".join(line.split("\t")[:3])
                if key in seen:
                    continue
                seen.add(key)
                parts = line.split("\t")
                # Ensure at least 4 columns: chrom start end gene
                if len(parts) < 4:
                    parts.append(gene)
                f.write("\t".join(parts[:4]) + "\n")
                count += 1
    return count


def _build_clinvar_variant_bed(
    disease_name: str,
    output: Path,
    vcf_path: Path = CLINVAR_VCF,
) -> int:
    """Build a BED of ClinVar P/LP variants whose phenotype matches the disease."""
    output.parent.mkdir(parents=True, exist_ok=True)
    if not vcf_path.exists():
        logger.warning("ClinVar VCF not found at %s; skipping ClinVar variant BED", vcf_path)
        return 0

    # Use canonical English disease key for keyword lookup so Chinese queries
    # do not produce empty normalized keywords that match every record.
    canonical = resolve_builtin_disease_key(disease_name)
    lookup_name = canonical or disease_name
    keywords = disease_keywords(lookup_name)
    # Drop empty/whitespace-only normalized keywords to avoid false positives.
    keywords = [kw for kw in keywords if _normalize(kw)]
    # Pre-filter to P/LP records only to avoid parsing the entire ClinVar VCF.
    proc = subprocess.run(
        [
            "bcftools", "view", "-H",
            '-i', 'CLNSIG~"pathogenic" || CLNSIG~"likely_pathogenic"',
            str(vcf_path),
        ],
        capture_output=True, text=True, check=False,
    )
    seen: set[str] = set()
    count = 0
    with open(output, "w", encoding="utf-8") as f:
        for line in proc.stdout.strip().split("\n"):
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _id, ref, alt, qual, filt, info = parts[:8]

            # Phenotype match
            norm_info = re.sub(r"[^a-z0-9 ]+", " ", info.lower())
            matched = False
            for kw in keywords:
                if _normalize(kw) in norm_info:
                    matched = True
                    break
            if not matched:
                continue

            m = re.search(r"CLNSIG=([^;]+)", info)
            sig = m.group(1).lower() if m else "pathogenic"
            start = int(pos) - 1
            end = start + len(ref)
            for a in alt.split(","):
                key = f"{chrom}\t{start}\t{end}\t{ref}>{a}"
                if key in seen:
                    continue
                seen.add(key)
                f.write(f"{key}\t{sig}\n")
                count += 1
    return count


def _build_gwas_bed(snps: list[dict], output: Path) -> int:
    """Write a BED of GWAS lead SNPs (1-based pos -> 0-based BED)."""
    output.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with open(output, "w", encoding="utf-8") as f:
        for s in snps:
            chrom = s.get("chrom", "")
            pos = int(s.get("pos", 0))
            if not chrom or pos <= 0:
                continue
            start = pos - 1
            end = pos
            gene = s.get("gene", "")
            rsid = s.get("rsid", "")
            note = s.get("note", "")
            f.write(f"{chrom}\t{start}\t{end}\t{rsid}\t{gene}\t{note}\n")
            count += 1
    return count


def _build_regions_bed(regions: list[dict], output: Path) -> int:
    """Write a BED of annotated regulatory/GWAS regions."""
    output.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with open(output, "w", encoding="utf-8") as f:
        for r in regions:
            chrom = r.get("chrom", "")
            start = int(r.get("start", 0))
            end = int(r.get("end", 0))
            if not chrom or start >= end:
                continue
            gene = r.get("gene", "")
            rtype = r.get("type", "")
            note = r.get("note", "")
            f.write(f"{chrom}\t{start}\t{end}\t{gene}\t{rtype}\t{note}\n")
            count += 1
    return count


def _normalize(text: str) -> str:
    return re.sub(r"[^a-z0-9 ]+", " ", text.lower()).strip()


class DiseaseReference:
    """Manages the on-disk cache of disease-specific reference data."""

    def __init__(self, disease_name: str, cache_root: Path = DISEASE_REF_CACHE) -> None:
        self.disease_name = disease_name
        self.canonical_key = resolve_builtin_disease_key(disease_name)
        # Use canonical key for cache directory when available, otherwise normalized name
        cache_key = self.canonical_key or _normalize_disease_key(disease_name)
        self.key = cache_key
        self.root = cache_root / cache_key
        self.metadata_path = self.root / "metadata.json"
        self.structured_json = self.root / "structured.json"
        self.genes_bed = self.root / "genes.bed"
        self.clinvar_bed = self.root / "clinvar_variants.bed"
        self.gwas_bed = self.root / "gwas_snps.bed"
        self.regions_bed = self.root / "regions.bed"
        self.literature_json = self.root / "literature.json"
        self.core_genes_txt = self.root / "core_genes.txt"

    def exists(self) -> bool:
        return self.metadata_path.exists()

    def build(self, force: bool = False) -> "DiseaseReference":
        """Build or refresh the disease reference cache."""
        if self.exists() and not force:
            logger.info("Disease reference cache found: %s", self.root)
            return self

        logger.info("Building disease reference cache for '%s' at %s", self.disease_name, self.root)
        self.root.mkdir(parents=True, exist_ok=True)

        builtin_key = self.canonical_key
        builtin = DISEASE_BUILTIN_REFS.get(builtin_key) if builtin_key else None

        # New structured gene_set (preferred) or legacy core_genes + gwas_loci_genes
        gene_set = builtin.get("gene_set", []) if builtin else []
        if gene_set:
            core_genes = {g["gene"] for g in gene_set if g.get("tier") in ("mendelian_high", "mendelian_mod")}
            gwas_loci_genes = {g["gene"] for g in gene_set if g.get("tier") not in ("mendelian_high", "mendelian_mod")}
        else:
            core_genes = set(builtin.get("core_genes", [])) if builtin else set()
            gwas_loci_genes = set(builtin.get("gwas_loci_genes", [])) if builtin else set()

        gwas_snps = builtin.get("gwas_lead_snps", []) if builtin else []
        regions = builtin.get("regions", []) if builtin else []
        literature = builtin.get("key_literature", []) if builtin else []

        # Persist structured reference data for downstream scoring/reporting
        structured = {
            "gene_set": gene_set,
            "gwas_lead_snps": gwas_snps,
            "regions": regions,
        }
        with open(self.structured_json, "w", encoding="utf-8") as f:
            json.dump(structured, f, indent=2, ensure_ascii=False)

        # Persist core genes (Mendelian/high-effect)
        with open(self.core_genes_txt, "w", encoding="utf-8") as f:
            for g in sorted(core_genes):
                f.write(f"{g}\n")

        # Build BEDs (core + GWAS loci)
        all_ref_genes = sorted(core_genes | gwas_loci_genes)
        genes_count = _build_genes_bed(all_ref_genes, self.genes_bed)
        clinvar_count = _build_clinvar_variant_bed(self.disease_name, self.clinvar_bed)
        gwas_count = _build_gwas_bed(gwas_snps, self.gwas_bed)
        regions_count = _build_regions_bed(regions, self.regions_bed)

        # Persist literature
        with open(self.literature_json, "w", encoding="utf-8") as f:
            json.dump(literature, f, indent=2, ensure_ascii=False)

        metadata = {
            "disease_name": self.disease_name,
            "key": self.key,
            "created": datetime.now(timezone.utc).isoformat(),
            "source": "builtin" if builtin else "generated",
            "counts": {
                "core_genes": len(core_genes),
                "gwas_loci_genes": len(gwas_loci_genes),
                "genes_bed_records": genes_count,
                "clinvar_variants": clinvar_count,
                "gwas_snps": gwas_count,
                "regions": regions_count,
                "literature_entries": len(literature),
            },
        }
        with open(self.metadata_path, "w", encoding="utf-8") as f:
            json.dump(metadata, f, indent=2, ensure_ascii=False)

        logger.info(
            "Disease reference built: %d genes, %d ClinVar variants, %d GWAS SNPs, %d regions, %d literature entries",
            len(core_genes), clinvar_count, gwas_count, regions_count, len(literature),
        )
        return self

    def load(self) -> dict:
        """Load reference data into a dict."""
        self.build()
        metadata = json.loads(self.metadata_path.read_text(encoding="utf-8"))
        core_genes = []
        if self.core_genes_txt.exists():
            core_genes = [line.strip() for line in self.core_genes_txt.read_text(encoding="utf-8").splitlines() if line.strip()]

        builtin_key = self.canonical_key
        builtin = DISEASE_BUILTIN_REFS.get(builtin_key) if builtin_key else None

        # Load structured reference data (new schema)
        structured = {"gene_set": [], "gwas_lead_snps": [], "regions": []}
        if self.structured_json.exists():
            structured = json.loads(self.structured_json.read_text(encoding="utf-8"))

        # Backward compatibility: derive gene_set from legacy fields if missing
        if not structured.get("gene_set") and builtin:
            legacy_core = set(builtin.get("core_genes", []))
            legacy_gwas = set(builtin.get("gwas_loci_genes", []))
            structured["gene_set"] = [
                {"gene": g, "tier": "mendelian_mod", "contribution_score": 0.5,
                 "penetrance": "moderate", "penetrance_score": 0.3, "evidence": "mixed", "note": ""}
                for g in sorted(legacy_core)
            ] + [
                {"gene": g, "tier": "gwas", "contribution_score": 0.15,
                 "penetrance": "very_low", "penetrance_score": 0.05, "evidence": "gwas", "note": ""}
                for g in sorted(legacy_gwas - legacy_core)
            ]

        gene_set = structured.get("gene_set", [])
        gwas_loci_genes = sorted({g["gene"] for g in gene_set if g.get("tier") not in ("mendelian_high", "mendelian_mod")})

        # Helper maps for scoring
        gene_contribution_map = {g["gene"]: g.get("contribution_score", 0.0) for g in gene_set}
        gene_penetrance_map = {g["gene"]: g.get("penetrance", "") for g in gene_set}
        gene_tier_map = {g["gene"]: g.get("tier", "gwas") for g in gene_set}
        snp_contribution_map = {s["rsid"]: s.get("contribution_score", 0.0) for s in structured.get("gwas_lead_snps", [])}

        literature = []
        if self.literature_json.exists():
            literature = json.loads(self.literature_json.read_text(encoding="utf-8"))

        return {
            "metadata": metadata,
            "core_genes": core_genes,
            "gwas_loci_genes": gwas_loci_genes,
            "gene_set": gene_set,
            "gwas_lead_snps": structured.get("gwas_lead_snps", []),
            "regions": structured.get("regions", []),
            "genes_bed": str(self.genes_bed),
            "clinvar_bed": str(self.clinvar_bed),
            "gwas_bed": str(self.gwas_bed),
            "regions_bed": str(self.regions_bed),
            "literature": literature,
            "gene_contribution_map": gene_contribution_map,
            "gene_penetrance_map": gene_penetrance_map,
            "gene_tier_map": gene_tier_map,
            "snp_contribution_map": snp_contribution_map,
        }


def get_disease_reference(disease_name: str, refresh: bool = False) -> dict:
    """Convenience wrapper to get or build a disease reference."""
    ref = DiseaseReference(disease_name)
    if refresh:
        ref.build(force=True)
    return ref.load()
