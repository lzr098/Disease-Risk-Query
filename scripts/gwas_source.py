"""GWAS data source for disease risk queries.

Primary source: built-in lead SNPs from constants.DISEASE_BUILTIN_REFS.
Optional extension: query the NHGRI-EBI GWAS Catalog REST API for additional
lead SNPs when a built-in list is absent or --gwas-online is enabled.
"""

from __future__ import annotations

import json
import logging
import re
import shutil
import tempfile
from pathlib import Path
from typing import Optional

from constants import DISEASE_BUILTIN_REFS, GRCH38_FASTA, resolve_builtin_disease_key

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


def _detect_chrom_prefix(vcf_path: Path) -> str:
    """Detect whether the VCF uses 'chr' prefixed contig names."""
    import subprocess
    header = subprocess.run(
        ["bcftools", "view", "-h", str(vcf_path)],
        capture_output=True, text=True, check=False,
    ).stdout
    for line in header.split("\n"):
        if line.startswith("##contig=<ID="):
            c = line.split("##contig=<ID=")[1].split(",")[0]
            return "chr" if c.startswith("chr") else ""
    return ""


def _fetch_fasta_ref(chrom: str, pos: int) -> str:
    """Fetch the reference allele at a GRCh38 position using samtools faidx.

    Returns empty string if the lookup fails.
    """
    import subprocess
    # FASTA uses unprefixed contig names (1, 2, X, Y, MT)
    c = chrom.lstrip("chr")
    region = f"{c}:{pos}-{pos}"
    proc = subprocess.run(
        ["samtools", "faidx", str(GRCH38_FASTA), region],
        capture_output=True, text=True, check=False,
    )
    lines = [ln.strip() for ln in proc.stdout.split("\n") if ln.strip() and not ln.startswith(">")]
    return lines[0].upper() if lines else ""


def _vcf_genotype_at(
    vcf_path: Path,
    chrom: str,
    pos: int,
    ref: str = "",
    alt: str = "",
    effect_allele: str = "",
) -> Optional[dict]:
    """Query a genotyped VCF for a SNP and return genotype info.

    If ref/alt are provided, returns the exact SNP. If not, returns the first
    SNP at the position whose REF/ALT alleles include the effect_allele.
    Returns None if no matching record is found.

    Note: This function spawns a separate bcftools process. For querying many
    positions, use _batch_query_vcf_positions instead.
    """
    import subprocess
    prefix = _detect_chrom_prefix(vcf_path)
    c = f"{prefix}{chrom.lstrip('chr')}"
    region = f"{c}:{pos}-{pos}"
    proc = subprocess.run(
        ["bcftools", "view", "-H", str(vcf_path), "-r", region],
        capture_output=True, text=True, check=False,
    )
    effect_allele = (effect_allele or "").upper()
    for line in proc.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 10:
            continue
        vcf_ref = parts[3].upper()
        vcf_alts = [a.upper() for a in parts[4].split(",")]
        # Exact match mode
        if ref and alt:
            if vcf_ref != ref.upper():
                continue
            if alt.upper() not in vcf_alts:
                continue
        else:
            # Flexible mode: require a SNP and that effect_allele matches REF or ALT
            if effect_allele and effect_allele not in ([vcf_ref] + vcf_alts):
                continue
            if not effect_allele:
                # No allele guidance; accept first biallelic SNP at this position
                if len(vcf_ref) != 1 or len(vcf_alts) != 1 or len(vcf_alts[0]) != 1:
                    continue
        matched_alt = alt.upper() if alt else (vcf_alts[0] if vcf_alts else "")
        fmt = parts[8].split(":")
        sample = parts[9].split(":")
        gt = sample[fmt.index("GT")] if "GT" in fmt else "./."
        return {
            "chrom": parts[0],
            "pos": int(parts[1]),
            "ref": vcf_ref,
            "alt": matched_alt,
            "gt": gt,
            "filter": parts[6],
            "format": parts[8],
            "sample": parts[9],
        }
    return None


def _normalize_chrom(chrom: str) -> str:
    return str(chrom).lstrip("chr")


def _parse_vcf_record_line(line: str) -> Optional[dict]:
    """Parse a single bcftools -H output line into a genotype dict."""
    parts = line.split("\t")
    if len(parts) < 10:
        return None
    vcf_ref = parts[3].upper()
    vcf_alts = [a.upper() for a in parts[4].split(",")]
    fmt = parts[8].split(":")
    sample = parts[9].split(":")
    gt = sample[fmt.index("GT")] if "GT" in fmt else "./."
    return {
        "chrom": parts[0],
        "pos": int(parts[1]),
        "ref": vcf_ref,
        "alts": vcf_alts,
        "gt": gt,
        "filter": parts[6],
        "format": parts[8],
        "sample": parts[9],
    }


def _match_record(
    records: list[dict],
    ref: str = "",
    alt: str = "",
    effect_allele: str = "",
) -> Optional[dict]:
    """Find a record matching ref/alt or effect_allele constraints."""
    ref_u = ref.upper()
    alt_u = alt.upper()
    effect_u = (effect_allele or "").upper()
    for rec in records:
        if ref and alt:
            if rec["ref"] != ref_u:
                continue
            if alt_u not in rec["alts"]:
                continue
            return {
                **rec,
                "alt": alt_u,
            }
        else:
            if effect_u and effect_u not in ([rec["ref"]] + rec["alts"]):
                continue
            if not effect_u:
                if len(rec["ref"]) != 1 or len(rec["alts"]) != 1 or len(rec["alts"][0]) != 1:
                    continue
            matched_alt = alt_u if alt else (rec["alts"][0] if rec["alts"] else "")
            return {**rec, "alt": matched_alt}
    return None


def _batch_query_vcf_positions(
    vcf_path: Path,
    positions: list[tuple[str, int]],
) -> dict[str, list[dict]]:
    """Query many VCF positions in a single bcftools call.

    Returns a dict mapping normalized "chrom:pos" to a list of record dicts.
    """
    import subprocess
    import tempfile
    from pathlib import Path

    result: dict[str, list[dict]] = {}
    if not positions:
        return result

    prefix = _detect_chrom_prefix(vcf_path)
    tmp_dir = Path(tempfile.mkdtemp(prefix="drq_gwas_batch_"))
    try:
        bed_path = tmp_dir / "positions.bed"
        with open(bed_path, "w", encoding="utf-8") as f:
            for chrom, pos in positions:
                c = f"{prefix}{_normalize_chrom(chrom)}"
                f.write(f"{c}\t{pos - 1}\t{pos}\n")

        proc = subprocess.run(
            ["bcftools", "view", "-H", "-R", str(bed_path), str(vcf_path)],
            capture_output=True, text=True, check=False,
        )
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    for line in proc.stdout.strip().split("\n"):
        if not line:
            continue
        rec = _parse_vcf_record_line(line)
        if rec is None:
            continue
        key = f"{_normalize_chrom(rec['chrom'])}:{rec['pos']}"
        result.setdefault(key, []).append(rec)
    return result


def check_gwas_lead_snps(
    vcf_path: Path,
    disease_name: str,
) -> list[dict]:
    """Check whether known GWAS lead SNPs are present in the sample.

    Returns a list of SNP dicts with a 'sample_gt' key added.
    For genotyped VCFs, any lead SNP not present in the VCF is treated as
    REF/REF (0/0) per project convention, not as a missing genotype.

    This function batches all position queries into a single bcftools call
    for performance.
    """
    snps = get_gwas_lead_snps(disease_name)
    if not snps:
        return []

    # Batch query all valid positions in one bcftools call
    valid_snps = []
    positions = []
    for s in snps:
        chrom = s.get("chrom", "")
        pos = int(s.get("pos", 0))
        if chrom and pos > 0:
            valid_snps.append((s, chrom, pos))
            positions.append((chrom, pos))

    position_records = _batch_query_vcf_positions(vcf_path, positions)

    results = []
    for s, chrom, pos in valid_snps:
        ref = s.get("ref", "")
        alt = s.get("alt", "")
        effect_allele = (s.get("effect_allele") or "").upper()
        other_allele = (s.get("other_allele") or "").upper()

        key = f"{_normalize_chrom(chrom)}:{pos}"
        records = position_records.get(key, [])
        gt = None
        if records:
            if ref and alt:
                gt = _match_record(records, ref=ref, alt=alt)
            elif effect_allele:
                gt = _match_record(records, effect_allele=effect_allele)
            else:
                gt = _match_record(records)

        if gt:
            results.append({**s, "sample_gt": gt})
            continue

        # SNP absent from VCF: infer REF/REF using FASTA reference allele
        inferred_ref = ref or _fetch_fasta_ref(chrom, pos)
        if not inferred_ref:
            results.append({**s, "sample_gt": None, "note": "ref_lookup_failed"})
            continue
        # Infer ALT as the non-reference allele among effect_allele/other_allele.
        candidates = {a for a in (effect_allele, other_allele) if a and a != inferred_ref}
        inferred_alt = alt or (candidates.pop() if candidates else "")
        results.append({
            **s,
            "sample_gt": {
                "chrom": chrom,
                "pos": pos,
                "ref": inferred_ref,
                "alt": inferred_alt,
                "gt": "0/0",
                "filter": ".",
                "format": "GT",
                "sample": "0/0",
                "inferred_ref_ref": True,
            },
            "note": "inferred_ref_ref",
        })

    # Append invalid SNPs (missing position) at the end
    for s in snps:
        chrom = s.get("chrom", "")
        pos = int(s.get("pos", 0))
        if not chrom or pos <= 0:
            results.append({**s, "sample_gt": None, "note": "missing_position"})

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
