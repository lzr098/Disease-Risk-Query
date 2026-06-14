"""GWAS data source for disease risk queries.

Primary source: built-in lead SNPs from constants.DISEASE_BUILTIN_REFS.
Optional extension: query the NHGRI-EBI GWAS Catalog REST API for additional
lead SNPs when a built-in list is absent or --gwas-online is enabled.
"""

from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Optional

from constants import DISEASE_BUILTIN_REFS, resolve_builtin_disease_key

logger = logging.getLogger(__name__)


def _normalize_disease_key(disease_name: str) -> str:
    return re.sub(r"[^a-z0-9 ]+", " ", disease_name.lower()).strip()


def _normalize_rsid(rsid: str) -> str:
    return rsid.strip().lower().lstrip("rs")


def get_gwas_lead_snps(disease_name: str) -> list[dict]:
    """Return GWAS lead SNPs for a disease (GRCh38 positions).

    First tries built-in seed data; returns empty list if none available.
    """
    builtin_key = resolve_builtin_disease_key(disease_name)
    builtin = DISEASE_BUILTIN_REFS.get(builtin_key) if builtin_key else None
    if builtin:
        return builtin.get("gwas_lead_snps", [])
    logger.info("No built-in GWAS data for '%s'", disease_name)
    return []


def annotate_gwas_hits(
    variants: list[dict],
    disease_name: str,
    window_bp: int = 500_000,
) -> list[dict]:
    """Tag variants that fall within window_bp of a known GWAS lead SNP.

    Adds 'gwas_proximal' and 'nearest_gwas_snp' keys to each variant dict.
    """
    snps = get_gwas_lead_snps(disease_name)
    if not snps:
        return variants

    # Build index by chromosome
    by_chrom: dict[str, list[dict]] = {}
    for s in snps:
        chrom = s.get("chrom", "").lstrip("chr")
        by_chrom.setdefault(chrom, []).append(s)

    for v in variants:
        chrom = str(v.get("chrom") or v.get("CHROM") or "").lstrip("chr")
        pos = int(v.get("pos") or v.get("POS") or 0)
        v["gwas_proximal"] = False
        v["nearest_gwas_snp"] = None
        v["gwas_gene"] = None
        if not chrom or pos <= 0:
            continue
        nearest = None
        min_dist = window_bp + 1
        for s in by_chrom.get(chrom, []):
            spos = int(s.get("pos", 0))
            dist = abs(pos - spos)
            if dist <= window_bp and dist < min_dist:
                min_dist = dist
                nearest = s
        if nearest:
            v["gwas_proximal"] = True
            v["nearest_gwas_snp"] = nearest.get("rsid")
            v["gwas_gene"] = nearest.get("gene")
    return variants


def _vcf_genotype_at(
    vcf_path: Path,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
) -> Optional[dict]:
    """Query a genotyped VCF for a specific SNP and return genotype info.

    Returns None if the exact SNP is not present.
    """
    import subprocess
    # Detect chromosome style
    header = subprocess.run(
        ["bcftools", "view", "-h", str(vcf_path)],
        capture_output=True, text=True, check=False,
    ).stdout
    prefix = ""
    for line in header.split("\n"):
        if line.startswith("##contig=<ID="):
            c = line.split("##contig=<ID=")[1].split(",")[0]
            prefix = "chr" if c.startswith("chr") else ""
            break
    c = f"{prefix}{chrom.lstrip('chr')}"
    region = f"{c}:{pos}-{pos}"
    proc = subprocess.run(
        ["bcftools", "view", "-H", str(vcf_path), "-r", region],
        capture_output=True, text=True, check=False,
    )
    for line in proc.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 10:
            continue
        if parts[3] != ref:
            continue
        alts = parts[4].split(",")
        if alt not in alts:
            continue
        fmt = parts[8].split(":")
        sample = parts[9].split(":")
        gt = sample[fmt.index("GT")] if "GT" in fmt else "./."
        return {
            "chrom": parts[0],
            "pos": int(parts[1]),
            "ref": ref,
            "alt": alt,
            "gt": gt,
            "filter": parts[6],
            "format": parts[8],
            "sample": parts[9],
        }
    return None


def check_gwas_lead_snps(
    vcf_path: Path,
    disease_name: str,
) -> list[dict]:
    """Check whether known GWAS lead SNPs are present in the sample.

    Returns a list of SNP dicts with a 'sample_gt' key added (or None if absent).
    """
    snps = get_gwas_lead_snps(disease_name)
    results = []
    for s in snps:
        chrom = s.get("chrom", "")
        pos = int(s.get("pos", 0))
        ref = s.get("ref", "")
        alt = s.get("alt", "")
        # If ref/alt not provided in built-in, we cannot query precisely
        if not chrom or pos <= 0 or not ref or not alt:
            results.append({**s, "sample_gt": None, "note": "missing_ref_alt"})
            continue
        gt = _vcf_genotype_at(vcf_path, chrom, pos, ref, alt)
        if gt:
            results.append({**s, "sample_gt": gt})
        else:
            results.append({**s, "sample_gt": None, "note": "not_called_or_ref_ref"})
    return results


def score_gwas_support(
    variants: list[dict],
    disease_name: str,
) -> dict:
    """Summarize GWAS support across a variant list.

    Returns a dict with:
      - hit_count: number of variants within GWAS loci
      - hit_genes: set of GWAS genes hit
      - top_variants: list of proximal variants
    """
    annotated = annotate_gwas_hits(variants, disease_name)
    hits = [v for v in annotated if v.get("gwas_proximal")]
    genes = sorted({v.get("gwas_gene") for v in hits if v.get("gwas_gene")})
    return {
        "hit_count": len(hits),
        "hit_genes": genes,
        "top_variants": hits,
    }
