"""Genome build liftover utilities using pyliftover + samtools.

Design choices:
- CrossMap is present in default venv metadata but its compiled dependency
  (pyBigWig) fails to load on macOS due to code-signature/Team-ID mismatch,
  making CrossMap unusable here.
- ``pyliftover`` is a lightweight pure-Python library that reads UCSC chain
  files and maps single coordinates. We use it for coordinate liftover and then
  rewrite the VCF ourselves with ``vcfpy`` (pure Python) and ``samtools faidx``.
- Chain files are interval-overlap binary mappings; storing them as SQL would be
  inefficient and lossy. Instead, we cache UCSC chain files locally under
  ``~/.workbuddy/data/liftover`` and reuse them.
"""

from __future__ import annotations

import gzip
import logging
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import vcfpy
from pyliftover import LiftOver
from vcfpy import Substitution

from constants import (
    GRCH37_ALIASES,
    GRCH38_ALIASES,
    GRCH38_FASTA,
    LIFTOVER_DIR,
    LIFTOVER_URLS,
)

logger = logging.getLogger(__name__)


def _to_ucsc_db(build: str) -> str:
    """Map a build name to UCSC database name."""
    b = build.lower()
    if b in {"grch37", "hg19"}:
        return "hg19"
    if b in {"grch38", "hg38"}:
        return "hg38"
    raise ValueError(f"Unsupported genome build for liftover: {build}")


def _chain_path(from_db: str, to_db: str) -> Path:
    """Expected local chain file path for a pyliftover direction."""
    return LIFTOVER_DIR / f"{from_db}To{to_db.capitalize()}.over.chain.gz"


def ensure_chain_file(
    source_build: str,
    target_build: str,
    chain_file: Optional[Path] = None,
) -> Path:
    """Ensure the liftover chain file exists, downloading via pyliftover if needed."""
    LIFTOVER_DIR.mkdir(parents=True, exist_ok=True)

    if chain_file is not None:
        path = Path(chain_file).expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(f"Custom chain file not found: {path}")
        # If a custom chain is provided but not in our cache dir, copy it there
        # under the expected name so pyliftover can find it with use_web=False.
        from_db = _to_ucsc_db(source_build)
        to_db = _to_ucsc_db(target_build)
        expected = _chain_path(from_db, to_db)
        if path != expected:
            shutil.copy2(path, expected)
        return expected

    from_db = _to_ucsc_db(source_build)
    to_db = _to_ucsc_db(target_build)
    expected = _chain_path(from_db, to_db)

    if expected.exists():
        return expected

    logger.info("Downloading liftover chain %s -> %s into %s", from_db, to_db, LIFTOVER_DIR)
    try:
        # pyliftover will fetch the chain file into cache_dir when use_web=True.
        LiftOver(from_db, to_db, cache_dir=str(LIFTOVER_DIR), use_web=True)
    except Exception as exc:
        url = LIFTOVER_URLS.get((from_db, to_db))
        raise RuntimeError(
            f"Failed to download liftover chain {from_db}->{to_db}: {exc}. "
            f"You may manually download {url} to {expected}."
        ) from exc

    if not expected.exists():
        url = LIFTOVER_URLS.get((from_db, to_db))
        raise RuntimeError(
            f"LiftOver chain file was not created at {expected}. "
            f"Please download {url} manually."
        )
    return expected


def detect_genome_build(vcf_path: Path) -> Optional[str]:
    """Infer genome build from VCF header.

    Looks for assembly= or reference= tags. Falls back to chromosome length
    heuristic for chr1 if header is ambiguous.
    """
    build: Optional[str] = None
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    try:
        with opener(vcf_path, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                if not line.startswith("#"):
                    break
                low = line.lower()
                if "grch38" in low or "hg38" in low:
                    build = "GRCh38"
                    break
                if "grch37" in low or "hg19" in low:
                    build = "GRCh37"
                    break
    except Exception as exc:
        logger.warning("Could not read VCF header for build detection: %s", exc)

    if build is None:
        try:
            result = subprocess.run(
                ["bcftools", "view", "-h", str(vcf_path)],
                capture_output=True, text=True, check=False,
            )
            for line in result.stdout.splitlines():
                if line.startswith("##contig=<ID=1") or line.startswith("##contig=<ID=chr1"):
                    if "length=" in line:
                        length_str = line.split("length=")[1].split(",")[0].split(">")[0]
                        length = int(length_str)
                        if length > 248_900_000:  # GRCh38 chr1 = 248,956,422
                            build = "GRCh38"
                        else:
                            build = "GRCh37"
                        break
        except Exception as exc:
            logger.warning("Contig length fallback failed: %s", exc)

    return build


def _fasta_has_chr_prefix(fasta_path: Path) -> bool:
    """Check whether the FASTA index uses chr-prefixed contig names."""
    fa_idx = fasta_path.with_suffix(fasta_path.suffix + ".fai")
    if not fa_idx.exists():
        return False
    with open(fa_idx, "r", encoding="utf-8") as f:
        return any(line.startswith("chr") for line in f)


def _normalize_chrom_for_fasta(chrom: str, fasta_has_chr: bool) -> str:
    """Match chromosome naming style (chr1 vs 1) to the FASTA index."""
    has_chr = chrom.startswith("chr")
    if fasta_has_chr and not has_chr:
        return f"chr{chrom}"
    if not fasta_has_chr and has_chr:
        return chrom[3:]
    return chrom


def _query_fasta(
    chrom: str,
    pos: int,
    length: int,
    fasta_path: Path,
    cache: Dict[Tuple[str, int, int], str],
    fasta_has_chr: Optional[bool] = None,
) -> str:
    """Return reference sequence (uppercase) at chrom:pos-pos+length-1."""
    key = (chrom, pos, length)
    if key in cache:
        return cache[key]

    if fasta_has_chr is None:
        fasta_has_chr = _fasta_has_chr_prefix(fasta_path)
    chrom_query = _normalize_chrom_for_fasta(chrom, fasta_has_chr)

    end = pos + length - 1
    proc = subprocess.run(
        ["samtools", "faidx", str(fasta_path), f"{chrom_query}:{pos}-{end}"],
        capture_output=True, text=True, check=False,
    )
    seq = "".join(proc.stdout.splitlines()[1:]).upper()
    cache[key] = seq
    return seq


def _normalize_chrom(chrom: str) -> Tuple[str, bool]:
    """Return (pyliftover_query_chrom, input_had_chr_prefix)."""
    if chrom.startswith("chr"):
        return chrom, True
    return f"chr{chrom}", False


def _invert_gt(gt: str) -> str:
    """Swap 0<->1 in a diploid genotype string, preserving separators."""
    if not gt or gt in (".", "./.", ".|."):
        return gt

    def _swap(allele: str) -> str:
        if allele == "0":
            return "1"
        if allele == "1":
            return "0"
        return allele

    if "|" in gt:
        return "|".join(_swap(a) for a in gt.split("|"))
    if "/" in gt:
        return "/".join(_swap(a) for a in gt.split("/"))
    return _swap(gt)


def liftover_vcf(
    input_vcf: Path,
    output_vcf: Path,
    source_build: str,
    target_build: str = "GRCh38",
    chain_file: Optional[Path] = None,
    fasta_path: Optional[Path] = None,
) -> dict:
    """Lift over a VCF between genome builds using pyliftover.

    Supports GRCh37/hg19 <-> GRCh38/hg38. For GRCh38 input this function copies
    the input to ``output_vcf`` unchanged and returns stats with
    ``lifted=False``.

    REF alleles are re-validated against the target FASTA after coordinate
    remapping; SNPs where the reference base changed between builds are
    re-encoded (REF/ALT swapped and the genotype inverted accordingly).
    Indels with a reference mismatch are dropped with a warning.

    Returns a dict with:
      - lifted (bool)
      - input_records, output_records, dropped
      - output_path
      - validation (dict from validate_ref_alleles)
    """
    source_build = source_build.lower()
    target_build = target_build.lower()
    fasta_path = fasta_path or GRCH38_FASTA

    if source_build in {b.lower() for b in GRCH38_ALIASES}:
        logger.info("Input already GRCh38; no liftover needed.")
        if input_vcf.resolve() != output_vcf.resolve():
            shutil.copy2(input_vcf, output_vcf)
            idx = input_vcf.with_suffix(input_vcf.suffix + ".csi")
            if idx.exists():
                shutil.copy2(idx, output_vcf.with_suffix(output_vcf.suffix + ".csi"))
        validation = validate_ref_alleles(output_vcf, fasta_path)
        return {
            "lifted": False,
            "input_records": 0,
            "output_records": 0,
            "dropped": 0,
            "output_path": str(output_vcf),
            "validation": validation,
        }

    if source_build not in {b.lower() for b in GRCH37_ALIASES}:
        raise ValueError(f"Unsupported source genome build: {source_build}")
    if target_build not in {b.lower() for b in GRCH38_ALIASES}:
        raise ValueError(f"Unsupported target genome build: {target_build}")

    ensure_chain_file(source_build, target_build, chain_file)
    from_db = _to_ucsc_db(source_build)
    to_db = _to_ucsc_db(target_build)
    lo = LiftOver(from_db, to_db, cache_dir=str(LIFTOVER_DIR), use_web=False)

    tmp_dir = Path(tempfile.mkdtemp(prefix="drq_liftover_"))
    try:
        reader = vcfpy.Reader.from_path(str(input_vcf))
        header = reader.header.copy()

        tmp_vcf = tmp_dir / "lifted.vcf"
        writer = vcfpy.Writer.from_path(str(tmp_vcf), header)

        fa_cache: Dict[Tuple[str, int, int], str] = {}
        fasta_has_chr = _fasta_has_chr_prefix(fasta_path)
        input_records = 0
        output_records = 0
        dropped = 0

        for record in reader:
            input_records += 1
            chrom = str(record.CHROM)
            query_chrom, had_chr = _normalize_chrom(chrom)

            # pyliftover uses 0-based coordinates
            mappings = lo.convert_coordinate(query_chrom, record.POS - 1)
            if not mappings:
                logger.debug("No liftover mapping for %s:%s", chrom, record.POS)
                dropped += 1
                continue

            # Prefer a forward-strand mapping; skip inversions
            target = None
            for m in mappings:
                if m[2] == "+":
                    target = m
                    break
            if target is None:
                logger.debug("Skipping inversion mapping for %s:%s", chrom, record.POS)
                dropped += 1
                continue

            target_chrom_ucsc = target[0]
            target_pos = target[1] + 1  # convert back to 1-based VCF POS
            target_chrom = target_chrom_ucsc[3:] if not had_chr else target_chrom_ucsc

            old_ref = str(record.REF).upper()
            old_alt = str(record.ALT[0].value).upper() if record.ALT else ""
            ref_len = len(old_ref)

            target_ref = _query_fasta(target_chrom, target_pos, ref_len, fasta_path, fa_cache, fasta_has_chr)
            if not target_ref:
                logger.warning("Could not query target FASTA at %s:%s", target_chrom, target_pos)
                dropped += 1
                continue

            swapped = False
            if old_ref == target_ref:
                new_ref = old_ref
                new_alt = old_alt
            elif len(old_ref) == 1 and len(old_alt) == 1 and old_alt == target_ref:
                # Reference allele changed between builds: swap REF/ALT and invert GT.
                new_ref = old_alt
                new_alt = old_ref
                swapped = True
            else:
                logger.debug(
                    "REF mismatch after liftover %s:%s -> %s:%s "
                    "(source REF=%s, target REF=%s); dropping record.",
                    chrom, record.POS, target_chrom, target_pos, old_ref, target_ref,
                )
                dropped += 1
                continue

            record.CHROM = target_chrom
            record.POS = target_pos
            record.REF = new_ref
            if record.ALT:
                record.ALT = [Substitution(record.ALT[0].type, new_alt)]

            if swapped:
                for call in record.calls:
                    gt = call.data.get("GT")
                    if gt is not None:
                        call.data["GT"] = _invert_gt(str(gt))
                record.update_calls(record.calls)

            writer.write_record(record)
            output_records += 1

        writer.close()
        reader.close()

        if output_records == 0:
            raise RuntimeError(
                "Liftover produced no usable records. Check that the input build "
                "and chain file are correct."
            )

        # Sort, compress, and index
        gz_path = output_vcf if str(output_vcf).endswith(".gz") else Path(str(output_vcf) + ".gz")
        subprocess.run(
            ["bcftools", "sort", "-Oz", "-o", str(gz_path), str(tmp_vcf)],
            check=True, capture_output=True, text=True,
        )
        subprocess.run(["bcftools", "index", str(gz_path)], check=True, capture_output=True, text=True)

        validation = validate_ref_alleles(gz_path, fasta_path)
        return {
            "lifted": True,
            "input_records": input_records,
            "output_records": output_records,
            "dropped": dropped,
            "output_path": str(gz_path),
            "validation": validation,
        }
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def validate_ref_alleles(
    vcf_path: Path,
    fasta_path: Path,
    sample_size: int = 100,
    max_mismatch_rate: float = 0.05,
) -> dict:
    """Sample variant sites and check REF allele matches FASTA."""
    import random

    cmd = ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\n", str(vcf_path)]
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    sites = []
    for line in result.stdout.strip().splitlines():
        parts = line.split("\t")
        if len(parts) == 3:
            chrom, pos, ref = parts
            sites.append((chrom, int(pos), ref))

    if not sites:
        return {"total_checked": 0, "mismatches": 0, "mismatch_rate": 0.0, "pass": True}

    sample = random.sample(sites, min(sample_size, len(sites)))
    mismatches = 0

    fa_cache: Dict[Tuple[str, int, int], str] = {}
    fasta_has_chr = _fasta_has_chr_prefix(fasta_path)
    for chrom, pos, ref in sample:
        try:
            fa_ref = _query_fasta(chrom, pos, len(ref), fasta_path, fa_cache, fasta_has_chr)
            if fa_ref != ref.upper():
                mismatches += 1
        except Exception as exc:
            logger.warning("FASTA validation failed at %s:%s: %s", chrom, pos, exc)
            mismatches += 1

    total = len(sample)
    rate = mismatches / total if total else 0.0
    return {
        "total_checked": total,
        "mismatches": mismatches,
        "mismatch_rate": rate,
        "pass": rate <= max_mismatch_rate,
    }
