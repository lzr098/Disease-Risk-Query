"""GPA / VEP 115 compatibility shim.

Problem:
  ``gpa-genomic-phenotype`` v0.10.5 invokes Docker VEP with
  ``--af_gnomad_exome`` and ``--af_gnomad_genome``. These flags were removed in
  VEP 115 (``ensemblorg/ensembl-vep:latest``). The modern equivalents are
  ``--af_gnomade`` (exome) and ``--af_gnomadg`` (genome), which produce
  ``gnomADe_AF`` / ``gnomADg_AF`` CSQ fields instead of the legacy
  ``gnomAD_AF``.

Solution:
  1. Pre-annotate the filtered VCF ourselves with the VEP 115-compatible flags.
  2. Monkey-patch ``dgra_input_parsers.VCFParser._csq_to_variant`` so that when
     the legacy ``gnomAD_AF`` is absent, it falls back to ``gnomADe_AF`` and
     then ``gnomADg_AF``.
  3. Pass the pre-annotated VCF to GPA. GPA detects ``INFO/CSQ`` and skips its
     own (incompatible) annotator.

This keeps all interpretation logic inside GPA while making the pipeline work
with the locally installed ``ensemblorg/ensembl-vep:latest`` image.
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional

from constants import GRCH38_FASTA, VEP_CACHE, VEP_IMAGE

logger = logging.getLogger(__name__)

_PATCHED = False
_VEP_FORKS = int(os.environ.get("DRQ_VEP_FORKS", "0"))


def _cpu_load_1m() -> float:
    """Return the 1-minute system load average. macOS/Linux compatible."""
    try:
        return os.getloadavg()[0]
    except Exception:
        return 99.0  # assume busy


def _cpu_idle_for_parallel() -> bool:
    """Return True if the system has idle capacity for parallel VEP forks.

    Heuristic: if the 1-minute load average is below 1.0 per available core
    (meaning < 100% CPU utilisation), parallelise with N-1 forks.
    """
    cores = os.cpu_count() or 1
    if cores < 2:
        return False
    load = _cpu_load_1m()
    # idle if total load < 1.0 per core (e.g. < 8 on 8-core machine)
    idle = load < float(cores) * 0.8  # 80% threshold
    logger.info("CPU check: 1m load %.2f, cores=%d, idle=%s", load, cores, idle)
    return idle


def _run_vep_single(part_vcf: Path, part_out: Path,
                    vep_args: List[str], volumes: List[str]) -> Path:
    """Run VEP Docker on a single VCF partition. Returns the output path."""
    cmd = ["docker", "run", "--rm", *volumes, VEP_IMAGE, "vep", *vep_args,
           "--input_file", f"/data/part/{part_vcf.name}",
           "--output_file", f"/data/part_out/{part_out.name}"]
    logger.debug("VEP partition: docker run --rm %s → %s", part_vcf.name, part_out.name)
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"VEP partition {part_vcf.name} failed: {proc.stderr[:500]}")
    return part_out


def _run_vep_parallel(partitions: List[Path], part_outputs: List[Path],
                      vep_args: List[str], volumes: List[str],
                      max_workers: int = 4) -> List[Path]:
    """Run VEP Docker in parallel across VCF partitions."""
    results: List[Path] = [None] * len(partitions)
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = {}
        for i, (pv, po) in enumerate(zip(partitions, part_outputs)):
            futures[pool.submit(_run_vep_single, pv, po, vep_args, volumes)] = i
        for future in as_completed(futures):
            i = futures[future]
            try:
                results[i] = future.result()
            except Exception as exc:
                logger.error("VEP partition %d failed: %s", i, exc)
                raise
    return results


def _split_vcf_for_parallel(vcf_path: Path, work_dir: Path, forks: int = 4) -> List[Path]:
    """Split a VCF into `forks` roughly equal partitions by chromosome.

    Tries ``bcftools index --stats`` first, then falls back to standard chromosomes.
    """
    import subprocess
    chroms: List[str] = []
    # Try to get chromosome list from index stats
    try:
        proc = subprocess.run(
            ["bcftools", "index", "--stats", str(vcf_path)],
            capture_output=True, text=True, check=False,
        )
        for line in proc.stdout.strip().split("\n"):
            if line.startswith("#") or not line.strip():
                continue
            chrom = line.split("\t")[0].strip()
            if chrom:
                chroms.append(chrom)
    except Exception:
        pass

    if not chroms:
        # Fallback: parse VCF header for contigs
        try:
            proc = subprocess.run(
                ["bcftools", "view", "-h", str(vcf_path)],
                capture_output=True, text=True, check=False,
            )
            for line in proc.stdout.strip().split("\n"):
                if line.startswith("##contig=<ID="):
                    import re
                    m = re.search(r"ID=([^,>]+)", line)
                    if m:
                        chroms.append(m.group(1))
        except Exception:
            pass

    if not chroms:
        chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]

    # Group chromosomes into forks partitions
    parts: List[List[str]] = [[] for _ in range(forks)]
    for i, c in enumerate(chroms):
        parts[i % forks].append(c)

    part_paths: List[Path] = []
    for i, chrom_group in enumerate(parts):
        if not chrom_group:
            continue
        part_path = work_dir / f"vep_part_{i}.vcf.gz"
        chrom_pattern = ",".join(chrom_group)
        subprocess.run(
            ["bcftools", "view", "-r", chrom_pattern, "-Oz", "-o", str(part_path), str(vcf_path)],
            capture_output=True, check=True,
        )
        subprocess.run(
            ["bcftools", "index", "-t", str(part_path)],
            capture_output=True, check=True,
        )
        part_paths.append(part_path)
        logger.info("VCF partition %d: %d chromosomes → %s", i, len(chrom_group), part_path.name)

    return part_paths


def patch_gpa_csq_parser() -> None:
    """Patch GPA's VCF parser to read gnomADe_AF / gnomADg_AF as gnomAD_AF."""
    global _PATCHED
    if _PATCHED:
        return

    try:
        from dgra_input_parsers import VCFParser
    except ImportError as exc:
        logger.warning("Could not import dgra_input_parsers for patching: %s", exc)
        return

    original = VCFParser._csq_to_variant

    def _compat(self, chrom: str, pos: int, ref: str, alt: str,
                csq: List[str], csq_map: Dict[str, int],
                gt: str, dp: str, gq: str, vaf: str) -> Dict[str, Any]:
        variant = original(self, chrom, pos, ref, alt, csq, csq_map, gt, dp, gq, vaf)

        if not variant.get("gnomAD_AF"):
            for field in ("gnomADe_AF", "gnomADg_AF"):
                idx = csq_map.get(field)
                if idx is not None and idx < len(csq):
                    val = csq[idx]
                    if val:
                        if "&" in val:
                            val = val.split("&")[0]
                        variant["gnomAD_AF"] = val
                        break
        return variant

    VCFParser._csq_to_variant = _compat
    _PATCHED = True
    logger.debug("Patched dgra_input_parsers.VCFParser for VEP 115 gnomAD fields")


def _vep_docker_available() -> bool:
    """Check whether the required VEP Docker image is present."""
    try:
        result = subprocess.run(
            ["docker", "image", "inspect", VEP_IMAGE],
            capture_output=True, text=True, check=False,
        )
        return result.returncode == 0
    except Exception:
        return False


def _vep_marker_path(output_vcf: Path) -> Path:
    """Return the path of the VEP checkpoint marker file."""
    return output_vcf.parent / f"{output_vcf.name}.vep_marker"


def _vep_output_fresh(input_vcf: Path, output_vcf: Path,
                      *, require_gnomad: bool = False) -> bool:
    """Check whether a valid VEP output exists that matches the current input.

    Returns True only when all conditions hold:
    1. The output VCF exists and is non-zero.
    2. The marker file exists and its stored fingerprint (mtime + size) matches
       the current input VCF.
    3. The output VCF passes ``gzip -t`` integrity check.
    4. If require_gnomad=True, the marker's gnomAD bit must be 1 (meaning the
       previous run included gnomAD AF flags).
    """
    if not output_vcf.exists() or output_vcf.stat().st_size == 0:
        return False

    marker = _vep_marker_path(output_vcf)
    if not marker.exists():
        return False

    try:
        parts = marker.read_text().strip().split()
        if len(parts) < 2:
            return False
        stored_mtime, stored_size = int(parts[0]), int(parts[1])
        gnomad_bit = parts[2] if len(parts) >= 3 else "0"
    except (ValueError, OSError):
        return False

    # Input fingerprint changed → stale output
    try:
        st = input_vcf.stat()
    except OSError:
        return False
    if stored_mtime != int(st.st_mtime) or stored_size != st.st_size:
        return False

    # gnomAD requirement check: if cache is now available but previous run
    # didn't include gnomAD AF, force re-run
    if require_gnomad and gnomad_bit != "1":
        return False

    # gzip integrity: detect truncated / corrupt output
    try:
        subprocess.run(
            ["gzip", "-t", str(output_vcf)],
            capture_output=True, check=True, text=False,
        )
    except (subprocess.CalledProcessError, OSError):
        return False

    return True


def _write_vep_marker(input_vcf: Path, output_vcf: Path,
                      *, gnomad_available: bool = False) -> None:
    """Persist the input fingerprint so the output can be reused later.

    Includes gnomAD availability status so that a run without gnomAD AF
    (cache unavailable) is automatically invalidated when the cache comes back.
    """
    try:
        st = input_vcf.stat()
        gnomad_bit = "1" if gnomad_available else "0"
        _vep_marker_path(output_vcf).write_text(
            f"{int(st.st_mtime)} {st.st_size} {gnomad_bit}"
        )
    except OSError:
        pass  # non-fatal — worst case VEP runs again next time


def pre_annotate_with_vep(
    input_vcf: Path,
    output_vcf: Path,
    genome: str = "GRCh38",
    vep_cache: Optional[Path] = None,
) -> Path:
    """Annotate a filtered VCF with VEP 115-compatible Docker flags.

    Uses the local VEP cache if available; falls back to database mode only
    when the cache is missing, and in that case omits gnomAD AF flags (which
    cannot be used with --database).
    """
    input_vcf = Path(input_vcf)
    output_vcf = Path(output_vcf)
    output_vcf.parent.mkdir(parents=True, exist_ok=True)

    cache_dir = Path(vep_cache or VEP_CACHE).expanduser()
    cache_available = cache_dir.exists() and (cache_dir / "homo_sapiens").exists()

    # --- Checkpoint resume: skip VEP if output is already fresh ---
    # require_gnomad=True → if previous run lacked gnomAD AF (cache was
    # unavailable at the time), force a fresh annotation so every report
    # variant carries a population frequency.
    if _vep_output_fresh(input_vcf, output_vcf, require_gnomad=cache_available):
        logger.info(
            "VEP checkpoint hit — output %s matches input %s, skipping annotation",
            output_vcf.name, input_vcf.name,
        )
        return output_vcf

    if not _vep_docker_available():
        raise RuntimeError(f"VEP Docker image {VEP_IMAGE} not found. Pull it with: docker pull {VEP_IMAGE}")

    assembly = genome if genome in {"GRCh38", "GRCh37"} else "GRCh38"

    # VEP writes uncompressed VCF; we compress afterwards.
    tmp_dir = Path(tempfile.mkdtemp(prefix="drq_vep_"))
    try:
        out_uncompressed = tmp_dir / "annotated.vcf"

        # Build core VEP flags (shared between single and parallel paths)
        core_flags = [
            "--vcf", "--assembly", assembly,
            "--canonical", "--mane", "--domains", "--protein", "--hgvs",
            "--numbers", "--check_existing", "--pubmed", "--symbol", "--biotype",
        ]

        if cache_available:
            core_flags.extend(["--cache", "--dir_cache", "/data/cache", "--offline"])
            core_flags.extend(["--af_gnomade", "--af_gnomadg"])
        else:
            core_flags.append("--database")
            if GRCH38_FASTA.exists():
                core_flags.extend(["--fasta", f"/data/fasta/{GRCH38_FASTA.name}"])

        use_parallel = _cpu_idle_for_parallel()
        forks = max(2, (os.cpu_count() or 4) - 1)

        if use_parallel:
            logger.info("CPU idle — running VEP with %d parallel forks", forks)
            part_dir = Path(tempfile.mkdtemp(prefix="drq_vep_parts_"))
            try:
                part_vcfs = _split_vcf_for_parallel(input_vcf, part_dir, forks)
                part_outputs = [part_dir / f"vep_out_{i}.vcf" for i in range(len(part_vcfs))]

                volumes_base = [
                    "-v", f"{cache_dir}:/data/cache:ro",
                    "-v", f"{part_dir}:/data/part:ro",
                    "-v", f"{part_dir}:/data/part_out",
                ]
                if not cache_available and GRCH38_FASTA.exists():
                    volumes_base.extend(["-v", f"{GRCH38_FASTA.parent}:/data/fasta:ro"])

                _run_vep_parallel(part_vcfs, part_outputs, core_flags, volumes_base, max_workers=forks)

                # Concatenate partition outputs (skip headers except first)
                with open(out_uncompressed, "w") as out_f:
                    for j, po in enumerate(part_outputs):
                        with open(po) as in_f:
                            for line in in_f:
                                if line.startswith("#"):
                                    if j == 0:
                                        out_f.write(line)
                                else:
                                    out_f.write(line)
                logger.info("Merged %d VEP partition outputs → %s", len(part_outputs), out_uncompressed)
            finally:
                shutil.rmtree(part_dir, ignore_errors=True)
        else:
            # Single-threaded VEP
            volumes = [
                "-v", f"{input_vcf.parent}:/data/input:ro",
                "-v", f"{tmp_dir}:/data/output",
            ]
            if cache_available:
                volumes.extend(["-v", f"{cache_dir}:/data/cache:ro"])
            elif not cache_available and GRCH38_FASTA.exists():
                volumes.extend(["-v", f"{GRCH38_FASTA.parent}:/data/fasta:ro"])

            cmd = ["docker", "run", "--rm", *volumes, VEP_IMAGE,
                   "vep",
                   "--input_file", f"/data/input/{input_vcf.name}",
                   "--output_file", f"/data/output/{out_uncompressed.name}",
                   *core_flags]
            logger.info("VEP 115: docker run --rm ... (%d variants)", _count_vcf_records(input_vcf))
            proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
            if proc.returncode != 0:
                raise RuntimeError(f"VEP 115 Docker failed: {proc.stderr[:1000]}")


        if not out_uncompressed.exists() or out_uncompressed.stat().st_size == 0:
            raise RuntimeError("VEP 115 produced no output file")

        # Compress and index
        gz_path = output_vcf if str(output_vcf).endswith(".gz") else Path(str(output_vcf) + ".gz")
        subprocess.run(
            ["bcftools", "view", "-Oz", "-o", str(gz_path), str(out_uncompressed)],
            check=True, capture_output=True, text=True,
        )
        subprocess.run(["bcftools", "index", str(gz_path)], check=True, capture_output=True, text=True)
        _write_vep_marker(input_vcf, gz_path, gnomad_available=cache_available)
        return gz_path
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def _count_vcf_records(vcf_path: Path) -> int:
    """Quick record count for logging."""
    try:
        result = subprocess.run(
            ["bcftools", "view", "-H", str(vcf_path)],
            capture_output=True, text=True, check=False,
        )
        return len(result.stdout.strip().splitlines())
    except Exception:
        return 0
