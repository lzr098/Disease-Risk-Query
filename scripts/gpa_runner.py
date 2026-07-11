"""Interface to gpa-genomic-phenotype skill for tiered variant classification.

This module enforces the GPA rule: do not use LLM knowledge; call the GPA
scripts for variant interpretation.
"""

from __future__ import annotations

import json
import logging
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Optional

from constants import DEFAULT_TISSUE, DISEASE_TISSUE_MAP
from gpa_compat import patch_gpa_csq_parser, pre_annotate_with_vep

logger = logging.getLogger(__name__)

GPA_SCRIPT = Path.home() / ".workbuddy" / "skills" / "dgra-genomic-risk" / "scripts" / "dgra_cli_wrapper.py"
PYTHON_BIN = Path.home() / ".workbuddy" / "binaries" / "python" / "envs" / "default" / "bin" / "python"


def infer_tissue(disease_name: str) -> str:
    """Map disease name to GPA tissue profile."""
    low = disease_name.lower()
    for tissue, keywords in DISEASE_TISSUE_MAP.items():
        if any(kw in low for kw in keywords):
            return tissue
    return DEFAULT_TISSUE


def _base_name(vcf: Path) -> str:
    name = vcf.name
    if name.endswith(".vcf.gz"):
        return name[:-7]
    if name.endswith(".vcf"):
        return name[:-4]
    return name


def _run_gpa_subprocess(
    annotated_vcf: Path,
    disease_name: str,
    disease_description: Optional[str],
    tissue: str,
    sex: str,
    age: Optional[int],
    offline: bool,
    spliceai: bool,
    two_phase: bool,
    progress_log: Optional[Path],
    output_json: Path,
) -> dict:
    """Run dgra_cli_wrapper.py in a subprocess to avoid asyncio event-loop conflicts.

    The dgra-genomic-risk skill uses asyncio.run() internally. Calling it from
    the same process as gpa-disease-risk-query can lead to 'coroutine never
    awaited' warnings and hard-to-debug hangs. Spawning a separate process
    eliminates this class of failures entirely.
    """
    if not GPA_SCRIPT.exists():
        raise FileNotFoundError(f"GPA wrapper not found: {GPA_SCRIPT}")
    if not PYTHON_BIN.exists():
        raise FileNotFoundError(f"Default Python interpreter not found: {PYTHON_BIN}")

    cmd = [
        str(PYTHON_BIN),
        str(GPA_SCRIPT),
        "--input-file", str(annotated_vcf),
        "--tissue", tissue,
        "--phenotypes", disease_name,
        "--disease-description", disease_description or disease_name,
        "--format", "vcf",
        "--output-json", str(output_json),
    ]
    if offline:
        cmd.append("--offline")
    if spliceai:
        cmd.append("--spliceai")
    if two_phase:
        cmd.append("--two-phase")
    if progress_log:
        cmd.extend(["--progress-log", str(progress_log)])

    logger.info("Running GPA subprocess: %s", " ".join(cmd))
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
            encoding="utf-8",
        )
    except Exception as exc:
        logger.exception("GPA subprocess failed to start")
        return {"success": False, "error": f"GPA subprocess failed to start: {exc}"}

    if result.returncode != 0:
        logger.error("GPA subprocess exited with code %s\nSTDERR: %s", result.returncode, result.stderr)
        return {
            "success": False,
            "error": f"GPA subprocess exited with code {result.returncode}: {result.stderr}",
        }

    try:
        with open(output_json, "r", encoding="utf-8") as f:
            gpa_result = json.load(f)
    except Exception as exc:
        logger.exception("Failed to read GPA subprocess output JSON")
        return {
            "success": False,
            "error": f"Failed to read GPA output JSON ({output_json}): {exc}",
            "stdout": result.stdout,
            "stderr": result.stderr,
        }

    return gpa_result


def run_gpa_on_filtered_vcf(
    filtered_vcf: Path,
    disease_name: str,
    disease_description: Optional[str] = None,
    tissue: Optional[str] = None,
    sex: str = "unknown",
    age: Optional[int] = None,
    offline: bool = False,
    spliceai: bool = True,
    two_phase: bool = False,
    progress_log: Optional[Path] = None,
) -> dict:
    """Run GPA analysis on a pre-filtered VCF.

    Uses a subprocess call to the dgra-genomic-risk CLI wrapper to avoid
    asyncio event-loop conflicts with the caller's process.

    Before GPA is invoked, the filtered VCF is pre-annotated with a VEP 115-
    compatible Docker command, and the GPA VCF parser is patched to read the
    modern ``gnomADe_AF`` / ``gnomADg_AF`` fields as ``gnomAD_AF``.
    """
    tissue = tissue or infer_tissue(disease_name)

    # Offline mode disables external API calls (SpliceAI, gnomAD REST, etc.)
    if offline and spliceai:
        logger.info("Offline mode: disabling SpliceAI external API queries")
        spliceai = False

    # Compatibility shim for VEP 115
    patch_gpa_csq_parser()
    annotated_vcf = filtered_vcf.parent / f"{_base_name(filtered_vcf)}.vep115.vcf.gz"
    annotated_vcf = pre_annotate_with_vep(filtered_vcf, annotated_vcf)

    output_json = filtered_vcf.parent / f"{_base_name(filtered_vcf)}.gpa_result.json"
    gpa_result = _run_gpa_subprocess(
        annotated_vcf=annotated_vcf,
        disease_name=disease_name,
        disease_description=disease_description,
        tissue=tissue,
        sex=sex,
        age=age,
        offline=offline,
        spliceai=spliceai,
        two_phase=two_phase,
        progress_log=progress_log,
        output_json=output_json,
    )

    if not gpa_result.get("success"):
        return gpa_result

    # Normalize summary fields
    results = gpa_result.get("results", {})
    summary = results.get("summary", {})
    tier1 = results.get("tier1_variants", [])
    tier2 = results.get("tier2_variants", [])
    tier3 = results.get("tier3_variants", [])
    multi_hit = results.get("multi_hit_details", results.get("multi_hit_genes", []))

    return {
        "success": True,
        "tissue": tissue,
        "tier1_variants": tier1,
        "tier2_variants": tier2,
        "tier3_variants": tier3,
        "multi_hit": multi_hit,
        "summary": summary,
        "report_md": gpa_result.get("report_md", ""),
        "raw": gpa_result,
    }


def save_gpa_json(result: dict, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, ensure_ascii=False, default=str)
    return path
