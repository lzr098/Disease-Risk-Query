"""Unified disease-space VCF query.

Builds a single BED covering all disease-associated genes and regions, queries
variants in that space, and separately queries known pathogenic/GWAS/PRS
variants with ref/ref inference.
"""

from __future__ import annotations

import logging
import math
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from disease_profile import DiseaseProfile, RegulatoryRegion, VariantWeight
from liftover import _query_fasta
from constants import GRCH38_FASTA

logger = logging.getLogger(__name__)


@dataclass
class KnownVariantGenotype:
    """Genotype result for a known variant."""

    variant: VariantWeight
    chrom: str
    pos: int
    ref: str
    alt: str
    gt: str
    dosage: int
    inferred_ref_ref: bool
    filter_status: str
    sample_format: str
    sample_values: str

    def to_dict(self):
        return {
            "variant": self.variant.to_dict(),
            "chrom": self.chrom,
            "pos": self.pos,
            "ref": self.ref,
            "alt": self.alt,
            "gt": self.gt,
            "dosage": self.dosage,
            "inferred_ref_ref": self.inferred_ref_ref,
            "filter_status": self.filter_status,
            "sample_format": self.sample_format,
            "sample_values": self.sample_values,
        }


def _count_vcf_records(vcf_path: Path) -> int:
    """Fast variant count using bcftools index -n when possible."""
    if not vcf_path.exists():
        return 0
    idx = vcf_path.with_suffix(vcf_path.suffix + ".csi")
    if not idx.exists():
        idx = vcf_path.with_suffix(vcf_path.suffix + ".tbi")
    if idx.exists():
        try:
            proc = subprocess.run(
                ["bcftools", "index", "-n", str(vcf_path)],
                capture_output=True, text=True, check=True,
            )
            return int(proc.stdout.strip())
        except Exception:
            pass
    proc = subprocess.run(
        ["bcftools", "view", "-H", str(vcf_path)],
        capture_output=True, text=True, check=False,
    )
    return proc.stdout.count("\n")


def _ensure_indexed(vcf_path: Path) -> None:
    csi = vcf_path.with_suffix(vcf_path.suffix + ".csi")
    tbi = vcf_path.with_suffix(vcf_path.suffix + ".tbi")
    if csi.exists() or tbi.exists():
        return
    logger.info("Indexing %s", vcf_path)
    subprocess.run(["bcftools", "index", str(vcf_path)], check=True)


def _build_unified_bed(
    profile: DiseaseProfile,
    gene_coords: dict[str, list[tuple[str, int, int, str]]],
    output: Path,
) -> int:
    """Write a unified BED of gene loci + regulatory regions.

    Returns the number of records written.
    """
    output.parent.mkdir(parents=True, exist_ok=True)
    seen: set[str] = set()
    count = 0
    with open(output, "w", encoding="utf-8") as f:
        # Gene loci
        for gene in profile.all_genes:
            intervals = gene_coords.get(gene, [])
            if not intervals:
                logger.warning("No coordinates for gene %s", gene)
                continue
            for chrom, start, end, strand in intervals:
                key = f"{chrom}\t{start}\t{end}"
                if key in seen:
                    continue
                seen.add(key)
                f.write(f"{chrom}\t{start}\t{end}\t{gene}\t0\t{strand}\n")
                count += 1

        # Regulatory / GWAS loci regions
        for r in profile.regulatory_regions:
            key = f"{r.chrom}\t{r.start}\t{r.end}"
            if key in seen:
                continue
            seen.add(key)
            gene = r.gene or "."
            f.write(f"{r.chrom}\t{r.start}\t{r.end}\t{gene}\t0\t.\t{r.rtype}\t{r.note}\n")
            count += 1

    logger.info("Wrote unified disease space BED: %d records -> %s", count, output)
    return count


def _bcftools_view_regions(vcf_path: Path, bed_path: Path, output: Path) -> None:
    """Run bcftools view -T to extract variants overlapping a BED."""
    output.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "bcftools", "view", "-T", str(bed_path),
        "-Oz", "-o", str(output),
        str(vcf_path),
    ]
    subprocess.run(cmd, check=True)
    _ensure_indexed(output)


def _parse_gt(fmt: str, sample: str, fmt_key: str = "GT") -> str:
    """Extract a format field from a sample column."""
    fmt_parts = fmt.split(":")
    sample_parts = sample.split(":")
    if fmt_key in fmt_parts and len(sample_parts) > fmt_parts.index(fmt_key):
        return sample_parts[fmt_parts.index(fmt_key)]
    return "./."


def _gt_to_dosage(gt: str, effect_allele: Optional[str], ref: str, alt: str) -> int:
    """Compute dosage of effect_allele from genotype string.

    Genotype alleles are VCF numeric indices (0=REF, 1=first ALT, etc.).
    """
    if not gt or gt in ("./.", ".|.", "."):
        return 0
    # Normalize phased/unphased separators
    indices = gt.replace("|", "/").split("/")
    if effect_allele:
        effect = effect_allele.upper()
        ref_upper = (ref or "").upper()
        alt_upper = (alt or "").upper()
        if effect == ref_upper:
            target = "0"
        elif effect == alt_upper:
            target = "1"
        else:
            # Effect allele not present at this site
            return 0
    else:
        # If no effect allele specified, count ALT copies
        target = "1"
    return sum(1 for a in indices if a == target)


def _normalize_allele(allele: str) -> str:
    return (allele or "").upper()


def _match_record(
    variant: VariantWeight,
    parts: list[str],
) -> Optional[KnownVariantGenotype]:
    """Match a VCF record to a VariantWeight and return genotype info."""
    chrom = parts[0]
    pos = int(parts[1])
    ref = parts[3].upper()
    alts = [a.upper() for a in parts[4].split(",")]

    # Allow position match even if ref/alt differ slightly (e.g. strand/normalization)
    if pos != variant.pos:
        return None

    # Prefer exact ref/alt match; fallback to effect allele match
    matched_alt = None
    if _normalize_allele(variant.ref) == ref and _normalize_allele(variant.alt) in alts:
        matched_alt = _normalize_allele(variant.alt)
    elif variant.effect_allele and variant.effect_allele.upper() in alts:
        matched_alt = variant.effect_allele.upper()
    elif variant.effect_allele and variant.effect_allele.upper() == ref:
        # Effect allele matches VCF REF — this is a REF/ALT swap relative
        # to the template.  Use the actual VCF ALT (the template's non-effect
        # allele) as matched_alt so the variant key is meaningful.
        matched_alt = alts[0] if alts else ref
    elif alts:
        matched_alt = alts[0]
    else:
        return None

    fmt = parts[8]
    sample = parts[9] if len(parts) > 9 else "."
    gt = _parse_gt(fmt, sample)
    dosage = _gt_to_dosage(gt, variant.effect_allele, ref, matched_alt)

    return KnownVariantGenotype(
        variant=variant,
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=matched_alt,
        gt=gt,
        dosage=dosage,
        inferred_ref_ref=False,
        filter_status=parts[6] if len(parts) > 6 else ".",
        sample_format=fmt,
        sample_values=sample,
    )


def _infer_ref_ref(variant: VariantWeight) -> KnownVariantGenotype:
    """Return a REF/REF genotype for a variant absent from the VCF.

    Dosage is computed from the VariantWeight structure:
    - REF = variant.ref (typically other_allele / non-effect allele)
    - ALT = variant.alt (typically effect_allele)
    - GT 0/0 = homozygous REF

    If effect_allele == REF, then 0/0 = 2 copies of the effect allele → dosage=2.
    If effect_allele == ALT, then 0/0 = 0 copies → dosage=0.
    If effect_allele matches neither, dosage=0 (unrecognised allele).
    """
    ref_allele = _normalize_allele(variant.ref)
    alt_allele = _normalize_allele(variant.alt)
    effect = (variant.effect_allele or "").upper()

    # Compute dosage: count of effect_allele in a 0/0 genotype
    if effect and effect == ref_allele:
        # Effect allele is the REF → 0/0 = homozygous for effect → 2 copies
        dosage = 2
    elif effect and effect == alt_allele:
        # Effect allele is the ALT → 0/0 = no ALT copies → 0 copies
        dosage = 0
    else:
        dosage = 0

    return KnownVariantGenotype(
        variant=variant,
        chrom=variant.chrom,
        pos=variant.pos,
        ref=ref_allele,
        alt=alt_allele,
        gt="0/0",
        dosage=dosage,
        inferred_ref_ref=True,
        filter_status=".",
        sample_format="GT",
        sample_values="0/0",
    )


def query_known_variants(
    vcf_path: Path,
    variants: list[VariantWeight],
) -> list[KnownVariantGenotype]:
    """Query exact positions of known variants, inferring ref/ref for absent ones."""
    if not variants:
        return []

    _ensure_indexed(vcf_path)

    # Build exact-position BED
    bed_lines: list[str] = []
    for v in variants:
        chrom = v.chrom.lstrip("chr") if not str(vcf_path).endswith(".vcf") else v.chrom
        # We keep the chrom as-is; bcftools will match whatever is in the VCF.
        bed_lines.append(f"{v.chrom}\t{v.pos - 1}\t{v.pos}\t{v.vcf_key}")

    with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as bed:
        bed.write("\n".join(bed_lines) + "\n")
        bed_path = Path(bed.name)

    try:
        proc = subprocess.run(
            ["bcftools", "view", "-H", "-T", str(bed_path), str(vcf_path)],
            capture_output=True, text=True, check=False,
        )
    finally:
        bed_path.unlink(missing_ok=True)

    found_keys: set[str] = set()
    results: list[KnownVariantGenotype] = []
    for line in proc.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 8:
            continue
        pos = int(parts[1])
        # Find matching VariantWeight
        for v in variants:
            if v.pos == pos:
                matched = _match_record(v, parts)
                if matched:
                    results.append(matched)
                    found_keys.add(v.vcf_key)
                break

    # Infer ref/ref for absent variants
    absent = [v for v in variants if v.vcf_key not in found_keys]
    for v in absent:
        results.append(_infer_ref_ref(v))

    logger.info(
        "Known variant query: %d found, %d inferred ref/ref",
        len(found_keys), len(absent),
    )
    return results


def query_disease_space(
    vcf_path: Path,
    profile: DiseaseProfile,
    gene_coords: dict[str, list[tuple[str, int, int, str]]],
    work_dir: Path,
) -> dict:
    """Query VCF within the unified disease space and known variants panel.

    Returns:
        dict with unified_bed, disease_space_vcf, known_variant_genotypes,
        and variant_count.
    """
    work_dir.mkdir(parents=True, exist_ok=True)

    unified_bed = work_dir / "unified_disease_space.bed"
    _build_unified_bed(profile, gene_coords, unified_bed)

    disease_space_vcf = work_dir / "disease_space_variants.vcf.gz"
    _bcftools_view_regions(vcf_path, unified_bed, disease_space_vcf)

    known_genotypes = query_known_variants(vcf_path, profile.all_known_variants)

    return {
        "unified_bed": unified_bed,
        "disease_space_vcf": disease_space_vcf,
        "known_variant_genotypes": known_genotypes,
        "variant_count": _count_vcf_records(disease_space_vcf),
    }
