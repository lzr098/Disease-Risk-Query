"""VCF filtering to disease gene regions plus ClinVar safety net."""

from __future__ import annotations

import gzip
import json
import logging
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Optional

from constants import (
    CLINVAR_PATHOGENIC_BED,
    GENCODE_GENE_LOCI_BED,
    GENCODE_GTF,
    GENE_COORDS_CACHE,
    VCF_COMMON_SNP_PRESENCE_THRESHOLD,
)
from clinvar_phenotype_matcher import build_disease_clinvar_bed
from gwas_source import get_gwas_lead_snps

logger = logging.getLogger(__name__)


def _detect_vcf_chrom_style(vcf_path: Path) -> str:
    """Return 'chr' if VCF contigs start with 'chr', else 'no_chr'."""
    result = subprocess.run(
        ["bcftools", "view", "-h", str(vcf_path)],
        capture_output=True, text=True, check=False,
    )
    for line in result.stdout.split("\n"):
        if line.startswith("##contig=<ID="):
            chrom = line.split("##contig=<ID=")[1].split(",")[0]
            return "chr" if chrom.startswith("chr") else "no_chr"
    return "no_chr"


def _normalize_bed_chrom(bed_path: Path, target_style: str, output_path: Path) -> Path:
    """Add or remove 'chr' prefix to match VCF style."""
    with open(bed_path, "r") as infile, open(output_path, "w") as outfile:
        for line in infile:
            stripped = line.strip()
            if not stripped or stripped.startswith("#") or stripped.startswith("track"):
                outfile.write(line)
                continue
            parts = stripped.split("\t")
            has_chr = parts[0].startswith("chr")
            if target_style == "chr" and not has_chr:
                parts[0] = f"chr{parts[0]}"
            elif target_style == "no_chr" and has_chr:
                parts[0] = parts[0][3:]
            outfile.write("\t".join(parts) + "\n")
    return output_path


def check_vcf_completeness(
    vcf_path: Path,
    disease_name: Optional[str] = None,
) -> dict:
    """Check whether common variants appear to have been filtered out.

    Uses disease-specific GWAS lead SNPs as anchor positions. A high-quality
    genotyping VCF should contain most of these positions even if the sample
    is REF/REF. If a large fraction are absent, the VCF was likely hard-filtered
    and the GWAS/common-variant dimension should be interpreted cautiously.

    Returns dict with presence rate, filtering flag, and total variant count.
    """
    anchors = get_gwas_lead_snps(disease_name or "") if disease_name else []
    if not anchors:
        return {
            "checked": False,
            "common_variants_filtered": False,
            "note": "No built-in GWAS anchor SNPs for this disease; cannot assess VCF completeness.",
        }

    tmp_dir = Path(tempfile.mkdtemp(prefix="drq_qc_"))
    try:
        # Build anchor BED (0-based half-open)
        bed_path = tmp_dir / "anchors.bed"
        with open(bed_path, "w", encoding="utf-8") as f:
            for s in anchors:
                chrom = s.get("chrom", "").lstrip("chr")
                pos = int(s.get("pos", 0))
                if chrom and pos > 0:
                    f.write(f"{chrom}\t{pos - 1}\t{pos}\t{s.get('rsid', '')}\n")

        vcf_style = _detect_vcf_chrom_style(vcf_path)
        norm_bed = tmp_dir / "anchors_norm.bed"
        _normalize_bed_chrom(bed_path, vcf_style, norm_bed)

        proc = subprocess.run(
            ["bcftools", "view", "-R", str(norm_bed), "-H", str(vcf_path)],
            capture_output=True, text=True, check=False,
        )

        hit_positions: set[str] = set()
        for line in proc.stdout.strip().split("\n"):
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            chrom = parts[0].lstrip("chr")
            hit_positions.add(f"{chrom}:{parts[1]}")

        expected_positions: set[str] = set()
        for s in anchors:
            chrom = s.get("chrom", "").lstrip("chr")
            pos = int(s.get("pos", 0))
            if chrom and pos > 0:
                expected_positions.add(f"{chrom}:{pos}")

        total = len(expected_positions)
        present = len(expected_positions & hit_positions)
        presence_rate = present / total if total else 0.0
        is_filtered = presence_rate < VCF_COMMON_SNP_PRESENCE_THRESHOLD

        # Total variant count: prefer index, fall back to line count
        total_variants = None
        idx_proc = subprocess.run(
            ["bcftools", "index", "-n", str(vcf_path)],
            capture_output=True, text=True, check=False,
        )
        if idx_proc.returncode == 0:
            try:
                total_variants = int(idx_proc.stdout.strip())
            except ValueError:
                total_variants = None
        if total_variants is None:
            count_proc = subprocess.run(
                ["bcftools", "view", "-H", str(vcf_path)],
                capture_output=True, text=True, check=False,
            )
            total_variants = count_proc.stdout.count("\n")

        return {
            "checked": True,
            "anchor_snps_checked": total,
            "anchor_snps_present": present,
            "presence_rate": round(presence_rate, 2),
            "common_variants_filtered": is_filtered,
            "total_variants": total_variants,
            "threshold": VCF_COMMON_SNP_PRESENCE_THRESHOLD,
            "note": (
                f"{present}/{total} anchor positions present ({presence_rate:.0%}); "
                f"common variants filtered; missing anchors treated as REF/REF"
                if is_filtered else "pass"
            ),
        }
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def build_gene_coordinates_cache(
    gtf_path: Path = GENCODE_GTF,
    cache_path: Path = GENE_COORDS_CACHE,
) -> dict[str, list[tuple[str, int, int, str]]]:
    """Parse GENCODE GTF and build gene -> coordinate intervals cache.

    Returns dict mapping gene symbol to list of (chrom, start, end, strand).
    Coordinates are 0-based half-open BED style.
    """
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    if cache_path.exists():
        try:
            with open(cache_path, "r", encoding="utf-8") as f:
                data = json.load(f)
            logger.info("Loaded gene coordinate cache with %d genes", len(data))
            return data
        except Exception as exc:
            logger.warning("Failed to load gene coords cache: %s", exc)

    logger.info("Building gene coordinate cache from %s (this may take ~30s)", gtf_path)
    gene_intervals: dict[str, list[tuple[str, int, int, str]]] = defaultdict(list)

    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    with opener(gtf_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = parts
            if feature != "gene":
                continue
            # Parse gene_name
            m = __import__("re").search(r'gene_name "([^"]+)"', attrs)
            if not m:
                continue
            gene = m.group(1)
            # GTF is 1-based inclusive
            gene_intervals[gene].append((chrom, int(start) - 1, int(end), strand))

    # Merge intervals per gene
    merged: dict[str, list[tuple[str, int, int, str]]] = {}
    for gene, intervals in gene_intervals.items():
        by_chrom: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
        for chrom, s, e, strand in intervals:
            by_chrom[chrom].append((s, e, strand))
        merged[gene] = []
        for chrom, ivs in by_chrom.items():
            ivs.sort()
            cur_start, cur_end, cur_strand = ivs[0]
            for s, e, strand in ivs[1:]:
                if s <= cur_end:
                    cur_end = max(cur_end, e)
                else:
                    merged[gene].append((chrom, cur_start, cur_end, cur_strand))
                    cur_start, cur_end, cur_strand = s, e, strand
            merged[gene].append((chrom, cur_start, cur_end, cur_strand))

    # JSON serialization: convert tuples to lists
    serializable = {
        gene: [[chrom, s, e, strand] for chrom, s, e, strand in intervals]
        for gene, intervals in merged.items()
    }
    with open(cache_path, "w", encoding="utf-8") as f:
        json.dump(serializable, f, indent=0)
    logger.info("Gene coordinate cache built: %d genes", len(serializable))
    return serializable


def build_gene_bed(
    genes: list[str],
    output_bed: Path,
    gtf_path: Path = GENCODE_GTF,
    coords_cache: Optional[dict] = None,
) -> Path:
    """Create a BED file covering all intervals for the requested genes."""
    cache = coords_cache if coords_cache is not None else build_gene_coordinates_cache(gtf_path)
    output_bed.parent.mkdir(parents=True, exist_ok=True)
    with open(output_bed, "w", encoding="utf-8") as f:
        written_genes = set()
        for gene in genes:
            intervals = cache.get(gene)
            if not intervals:
                continue
            written_genes.add(gene)
            for interval in intervals:
                chrom, start, end, strand = interval
                f.write(f"{chrom}\t{start}\t{end}\t{gene}\t0\t{strand}\n")
    logger.info("Wrote BED for %d/%d genes: %s", len(written_genes), len(genes), output_bed)
    return output_bed


def filter_vcf_by_gene_set(
    input_vcf: Path,
    output_vcf: Path,
    genes: list[str],
    clinvar_safetynet: bool = True,
    gtf_path: Path = GENCODE_GTF,
    disease_name: Optional[str] = None,
) -> dict:
    """Filter VCF to disease gene regions + optional ClinVar P/LP safety net.

    If disease_name is provided, the ClinVar safety net is filtered to keep only
    variants whose phenotype matches the disease, reducing unrelated pathogenic
    hits (e.g. FLG ichthyosis variants in an Alzheimer query).

    Returns dict with paths and counts.
    """
    # Load gene coordinates once and reuse for BED building and stats
    coords_cache = build_gene_coordinates_cache(gtf_path)

    tmp_dir = Path(tempfile.mkdtemp(prefix="drq_filter_"))
    try:
        # Build gene BED
        gene_bed = tmp_dir / "genes.bed"
        build_gene_bed(genes, gene_bed, gtf_path, coords_cache=coords_cache)

        # Normalize chromosome style to match VCF
        vcf_style = _detect_vcf_chrom_style(input_vcf)
        norm_gene_bed = tmp_dir / "genes_norm.bed"
        _normalize_bed_chrom(gene_bed, vcf_style, norm_gene_bed)

        # Sort BED
        sorted_gene_bed = tmp_dir / "genes_sorted.bed"
        subprocess.run(
            ["sort", "-k1,1", "-k2,2n", "-o", str(sorted_gene_bed), str(norm_gene_bed)],
            check=True, capture_output=True, text=True,
        )

        # Region filter
        region_vcf = tmp_dir / "region.vcf.gz"
        subprocess.run(
            ["bcftools", "view", "-T", str(sorted_gene_bed), "-Oz", "-o", str(region_vcf), str(input_vcf)],
            check=True, capture_output=True, text=True,
        )

        vcfs_to_merge = [region_vcf]

        # ClinVar safety net (disease-matched when disease_name is given)
        # When a disease is specified, we only use ClinVar variants whose
        # phenotype names match that disease. Falling back to the broad
        # all-pathogenic BED pulls in P/LP variants for unrelated disorders
        # (e.g. VWF, CD36, MYH11) and produces false-positive Tier 1/2 hits.
        clinvar_count = 0
        if clinvar_safetynet:
            if disease_name:
                disease_clinvar_bed = tmp_dir / "clinvar_disease_matched.bed"
                try:
                    build_disease_clinvar_bed(disease_name, disease_clinvar_bed)
                    if disease_clinvar_bed.exists() and disease_clinvar_bed.stat().st_size > 0:
                        source_bed = disease_clinvar_bed
                    else:
                        logger.warning(
                            "Disease-matched ClinVar BED is empty for '%s'; "
                            "skipping ClinVar safety net to avoid unrelated "
                            "pathogenic variants.",
                            disease_name,
                        )
                        source_bed = None
                except Exception as exc:
                    logger.warning(
                        "Failed to build disease-matched ClinVar BED for '%s'; "
                        "skipping ClinVar safety net: %s",
                        disease_name, exc,
                    )
                    source_bed = None
            elif CLINVAR_PATHOGENIC_BED.exists():
                source_bed = CLINVAR_PATHOGENIC_BED
            else:
                source_bed = None

            if source_bed is not None:
                norm_clinvar_bed = tmp_dir / "clinvar_norm.bed"
                _normalize_bed_chrom(source_bed, vcf_style, norm_clinvar_bed)
                clinvar_vcf = tmp_dir / "clinvar.vcf.gz"
                subprocess.run(
                    ["bcftools", "view", "-T", str(norm_clinvar_bed), "-Oz", "-o", str(clinvar_vcf), str(input_vcf)],
                    check=True, capture_output=True, text=True,
                )
                clinvar_count = int(subprocess.run(
                    ["bcftools", "view", "-H", str(clinvar_vcf)],
                    capture_output=True, text=True, check=False,
                ).stdout.count("\n"))
                if clinvar_count > 0:
                    vcfs_to_merge.append(clinvar_vcf)

        output_vcf.parent.mkdir(parents=True, exist_ok=True)
        if len(vcfs_to_merge) == 1:
            shutil.copy2(region_vcf, output_vcf)
        else:
            # Merge and deduplicate
            for vcf in vcfs_to_merge:
                subprocess.run(["bcftools", "index", str(vcf)], check=True, capture_output=True, text=True)
            concat = tmp_dir / "concat.vcf.gz"
            subprocess.run(
                ["bcftools", "concat", "-a"] + [str(v) for v in vcfs_to_merge] + ["-Oz", "-o", str(concat)],
                check=True, capture_output=True, text=True,
            )
            sorted_merged = tmp_dir / "sorted.vcf.gz"
            subprocess.run(
                ["bcftools", "sort", str(concat), "-Oz", "-o", str(sorted_merged)],
                check=True, capture_output=True, text=True,
            )
            subprocess.run(
                ["bcftools", "norm", "-d", "snps", str(sorted_merged), "-Oz", "-o", str(output_vcf)],
                check=True, capture_output=True, text=True,
            )

        # Index output
        subprocess.run(["bcftools", "index", str(output_vcf)], check=True, capture_output=True, text=True)

        region_count = int(subprocess.run(
            ["bcftools", "view", "-H", str(region_vcf)],
            capture_output=True, text=True, check=False,
        ).stdout.count("\n"))
        total_count = int(subprocess.run(
            ["bcftools", "view", "-H", str(output_vcf)],
            capture_output=True, text=True, check=False,
        ).stdout.count("\n"))

        return {
            "output_vcf": output_vcf,
            "region_variants": region_count,
            "clinvar_variants": clinvar_count,
            "total_variants": total_count,
            "genes_with_coords": len(set(genes) & set(coords_cache.keys())),
        }
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)
