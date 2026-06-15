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
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional

from constants import GRCH38_FASTA, VEP_CACHE, VEP_IMAGE

logger = logging.getLogger(__name__)

_PATCHED = False


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


def _vep_output_fresh(input_vcf: Path, output_vcf: Path) -> bool:
    """Check whether a valid VEP output exists that matches the current input.

    Returns True only when all three conditions hold:
    1. The output VCF exists and is non-zero.
    2. The marker file exists and its stored fingerprint (mtime + size) matches
       the current input VCF.
    3. The output VCF passes ``gzip -t`` integrity check.
    """
    if not output_vcf.exists() or output_vcf.stat().st_size == 0:
        return False

    marker = _vep_marker_path(output_vcf)
    if not marker.exists():
        return False

    try:
        parts = marker.read_text().strip().split()
        if len(parts) != 2:
            return False
        stored_mtime, stored_size = int(parts[0]), int(parts[1])
    except (ValueError, OSError):
        return False

    # Input fingerprint changed → stale output
    try:
        st = input_vcf.stat()
    except OSError:
        return False
    if stored_mtime != int(st.st_mtime) or stored_size != st.st_size:
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


def _write_vep_marker(input_vcf: Path, output_vcf: Path) -> None:
    """Persist the input fingerprint so the output can be reused later."""
    try:
        st = input_vcf.stat()
        _vep_marker_path(output_vcf).write_text(f"{int(st.st_mtime)} {st.st_size}")
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

    # --- Checkpoint resume: skip VEP if output is already fresh ---
    if _vep_output_fresh(input_vcf, output_vcf):
        logger.info(
            "VEP checkpoint hit — output %s matches input %s, skipping annotation",
            output_vcf.name, input_vcf.name,
        )
        return output_vcf

    if not _vep_docker_available():
        raise RuntimeError(f"VEP Docker image {VEP_IMAGE} not found. Pull it with: docker pull {VEP_IMAGE}")

    cache_dir = Path(vep_cache or VEP_CACHE).expanduser()
    cache_available = cache_dir.exists() and (cache_dir / "homo_sapiens").exists()

    assembly = genome if genome in {"GRCh38", "GRCh37"} else "GRCh38"

    # VEP writes uncompressed VCF; we compress afterwards.
    tmp_dir = Path(tempfile.mkdtemp(prefix="drq_vep_"))
    try:
        out_uncompressed = tmp_dir / "annotated.vcf"

        volumes = [
            "-v", f"{input_vcf.parent}:/data/input:ro",
            "-v", f"{tmp_dir}:/data/output",
        ]
        vep_args = [
            "vep",
            "--input_file", f"/data/input/{input_vcf.name}",
            "--output_file", f"/data/output/{out_uncompressed.name}",
            "--vcf",
            "--assembly", assembly,
            "--canonical",
            "--mane",
            "--domains",
            "--protein",
            "--hgvs",
            "--numbers",
            "--check_existing",
            "--pubmed",
            "--symbol",
            "--biotype",
        ]

        if cache_available:
            volumes.extend(["-v", f"{cache_dir}:/data/cache:ro"])
            vep_args.extend([
                "--cache",
                "--dir_cache", "/data/cache",
                "--offline",
                # VEP 115 gnomAD AF flags
                "--af_gnomade",
                "--af_gnomadg",
            ])
            logger.info("Running VEP 115 in offline/cache mode with gnomAD exome/genome AF")
        else:
            logger.warning(
                "VEP cache not found at %s; falling back to --database mode "
                "(gnomAD AF will be unavailable).", cache_dir
            )
            vep_args.append("--database")
            # Mount GRCh38 FASTA for database mode if available
            if GRCH38_FASTA.exists():
                volumes.extend(["-v", f"{GRCH38_FASTA.parent}:/data/fasta:ro"])
                vep_args.extend([
                    "--fasta", f"/data/fasta/{GRCH38_FASTA.name}",
                ])

        cmd = ["docker", "run", "--rm", *volumes, VEP_IMAGE, *vep_args]
        logger.info("VEP 115 command: docker run --rm ... (%d variants)",
                    _count_vcf_records(input_vcf))

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
        _write_vep_marker(input_vcf, gz_path)
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
