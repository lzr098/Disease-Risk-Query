"""APOE genotype check for Alzheimer disease risk queries.

For a genotyped VCF, absence of a variant at the APOE ε2/ε4 sites means the
sample is homozygous reference at those positions. This module reports the
inferred APOE genotype when possible, or flags when the sites are not callable.
"""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path
from typing import Optional

from constants import APOE_GRCH38

logger = logging.getLogger(__name__)


def _detect_chrom_style(vcf_path: Path) -> str:
    """Return 'chr' if VCF contigs use chr prefix, else ''."""
    result = subprocess.run(
        ["bcftools", "view", "-h", str(vcf_path)],
        capture_output=True, text=True, check=False,
    )
    for line in result.stdout.split("\n"):
        if line.startswith("##contig=<ID="):
            chrom = line.split("##contig=<ID=")[1].split(",")[0]
            return "chr" if chrom.startswith("chr") else ""
    return ""


def _region_string(chrom: str, start: int, end: int, chr_prefix: str) -> str:
    c = f"{chr_prefix}{chrom.lstrip('chr')}"
    return f"{c}:{start}-{end}"


def _vcf_records(vcf_path: Path, region: str) -> list[dict]:
    """Return parsed VCF records for a region."""
    result = subprocess.run(
        ["bcftools", "view", "-H", str(vcf_path), "-r", region],
        capture_output=True, text=True, check=False,
    )
    records = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 8:
            continue
        records.append({
            "chrom": parts[0],
            "pos": int(parts[1]),
            "id": parts[2],
            "ref": parts[3],
            "alt": parts[4],
            "qual": parts[5],
            "filter": parts[6],
            "info": parts[7],
            "format": parts[8] if len(parts) > 8 else None,
            "sample": parts[9] if len(parts) > 9 else None,
        })
    return records


def _allele_at(records: list[dict], pos: int, ref: str, alt: str) -> Optional[str]:
    """Return genotype allele for a specific SNP, or None if not present."""
    for rec in records:
        if rec["pos"] == pos and rec["ref"] == ref and alt in rec["alt"].split(","):
            if not rec["format"] or not rec["sample"]:
                return None
            fmt = rec["format"].split(":")
            sample = rec["sample"].split(":")
            gt_idx = fmt.index("GT") if "GT" in fmt else -1
            if gt_idx < 0 or gt_idx >= len(sample):
                return None
            gt = sample[gt_idx]
            # Normalize phased/unphased
            alleles = gt.replace("|", "/").split("/")
            alt_idx = rec["alt"].split(",").index(alt) + 1
            counts = sum(1 for a in alleles if a == str(alt_idx))
            if counts == 0:
                return "ref"
            if counts == 1:
                return "het"
            return "hom"
    return None


def check_apoe(
    vcf_path: Path,
    require_callable: bool = True,
) -> dict:
    """Infer APOE genotype from a GRCh38 genotyped VCF.

    Returns a dict with:
      - present: bool (whether either APOE site was found)
      - rs7412_status, rs429358_status: 'ref'/'het'/'hom'/None
      - inferred_allele: 'ε2/ε2', 'ε2/ε3', 'ε2/ε4', 'ε3/ε3', 'ε3/ε4', 'ε4/ε4',
        or None if sites are absent
      - warning: human-readable message if APOE sites are missing
      - is_reference: True if both sites are reference (ε3/ε3)
      - region: region checked
    """
    prefix = _detect_chrom_style(vcf_path)
    region = _region_string(
        APOE_GRCH38["chrom"], APOE_GRCH38["start"], APOE_GRCH38["end"], prefix
    )
    records = _vcf_records(vcf_path, region)

    rs7412 = APOE_GRCH38["rs7412"]
    rs429358 = APOE_GRCH38["rs429358"]

    status_7412 = _allele_at(records, rs7412["pos"], rs7412["ref"], rs7412["alt"])
    status_429358 = _allele_at(records, rs429358["pos"], rs429358["ref"], rs429358["alt"])

    present = status_7412 is not None or status_429358 is not None

    result = {
        "present": present,
        "rs7412_status": status_7412,
        "rs429358_status": status_429358,
        "inferred_allele": None,
        "is_reference": False,
        "warning": None,
        "region": region,
    }

    if not present:
        result["warning"] = (
            "APOE ε2/ε4 sites (rs7412/rs429358) were not found in the VCF. "
            "This VCF appears to have been filtered for rare variants; APOE "
            "genotype cannot be determined and late-onset Alzheimer risk "
            "assessment is incomplete."
        )
        return result

    # Infer alleles from the two SNPs.
    # ε4: rs429358-T + rs7412-T  (T at both)
    # ε3: rs429358-C + rs7412-T  (C at rs429358, T at rs7412)
    # ε2: rs429358-C + rs7412-C  (C at both)
    # We encode presence (het/hom) of the ALT allele at each site.
    has_t_7412 = status_7412 in ("ref",)  # reference T
    has_c_7412 = status_7412 in ("het", "hom")  # alternate C
    has_t_429358 = status_429358 in ("het", "hom")  # alternate T
    has_c_429358 = status_429358 in ("ref",)  # reference C

    def _alleles():
        # For each chromosome, determine ε allele based on the two positions.
        # We need phased GT ideally; for unphased we enumerate possible alleles.
        # Simplified: use counts of each site to infer diplotype.
        counts = {"ε2": 0, "ε3": 0, "ε4": 0}
        if status_7412 == "hom" and status_429358 == "hom":
            # Both sites homozygous -> unambiguous
            if has_c_7412 and has_t_429358:
                counts["ε4"] = 2
            elif has_t_7412 and has_c_429358:
                counts["ε3"] = 2
            elif has_c_7412 and has_c_429358:
                counts["ε2"] = 2
        else:
            # Count chromosomes carrying each allele combination
            # rs7412 C count
            c7412 = 0
            if status_7412 == "het":
                c7412 = 1
            elif status_7412 == "hom":
                c7412 = 2
            # rs429358 T count
            t429358 = 0
            if status_429358 == "het":
                t429358 = 1
            elif status_429358 == "hom":
                t429358 = 2
            # Each chromosome with T at rs429358 is ε4
            counts["ε4"] = t429358
            # Remaining chromosomes with C at rs7412 are ε2
            counts["ε2"] = max(0, c7412 - t429358)
            # Rest are ε3
            counts["ε3"] = 2 - counts["ε4"] - counts["ε2"]

        # Build diplotype string
        alleles = []
        for a in ("ε2", "ε3", "ε4"):
            alleles.extend([a] * counts[a])
        return "/".join(alleles) if len(alleles) == 2 else None

    inferred = _alleles()
    result["inferred_allele"] = inferred
    result["is_reference"] = inferred == "ε3/ε3"
    if result["is_reference"]:
        result["warning"] = (
            "APOE genotype inferred as ε3/ε3 (reference); this is the most "
            "common genotype and does not carry the ε4 risk allele."
        )
    return result


if __name__ == "__main__":
    import sys
    import json
    vcf = Path(sys.argv[1])
    print(json.dumps(check_apoe(vcf), indent=2, ensure_ascii=False))
