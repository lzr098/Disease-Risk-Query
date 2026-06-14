"""Interface to gpa-genomic-phenotype skill for tiered variant classification.

This module enforces the GPA rule: do not use LLM knowledge; call the GPA
scripts for variant interpretation.
"""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path
from typing import Optional

from constants import DEFAULT_TISSUE, DISEASE_TISSUE_MAP
from gpa_compat import patch_gpa_csq_parser, pre_annotate_with_vep

logger = logging.getLogger(__name__)

GPA_SCRIPT = Path.home() / ".workbuddy" / "skills" / "dgra-genomic-risk" / "scripts" / "dgra_cli_wrapper.py"


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


def run_gpa_on_filtered_vcf(
    filtered_vcf: Path,
    disease_name: str,
    tissue: Optional[str] = None,
    sex: str = "unknown",
    age: Optional[int] = None,
    offline: bool = False,
    spliceai: bool = True,
    two_phase: bool = False,
    progress_log: Optional[Path] = None,
) -> dict:
    """Run GPA analysis on a pre-filtered VCF.

    Uses the direct Python API path by importing dgra_cli_wrapper functions
    to avoid subprocess overhead and keep logs within the same process.

    Before GPA is invoked, the filtered VCF is pre-annotated with a VEP 115-
    compatible Docker command, and the GPA VCF parser is patched to read the
    modern ``gnomADe_AF`` / ``gnomADg_AF`` fields as ``gnomAD_AF``.
    """
    tissue = tissue or infer_tissue(disease_name)
    if not GPA_SCRIPT.exists():
        raise FileNotFoundError(f"GPA wrapper not found: {GPA_SCRIPT}")

    gpa_dir = GPA_SCRIPT.parent
    sys.path.insert(0, str(gpa_dir))

    try:
        import dgra_cli_wrapper as wrapper

        # Compatibility shim for VEP 115
        patch_gpa_csq_parser()
        annotated_vcf = filtered_vcf.parent / f"{_base_name(filtered_vcf)}.vep115.vcf.gz"
        annotated_vcf = pre_annotate_with_vep(filtered_vcf, annotated_vcf)

        result = wrapper.run_gpa_from_file(
            input_path=annotated_vcf,
            tissue=tissue,
            user_phenotypes=disease_name,
            offline=offline,
            somatic=False,
            spliceai_enabled=spliceai,
            disease_description=disease_name,
            annotator="auto",
            two_phase=two_phase,
            progress_log_path=str(progress_log) if progress_log else None,
        )
    except Exception as exc:
        logger.exception("GPA analysis failed")
        return {"success": False, "error": f"GPA analysis failed: {exc}"}
    finally:
        if str(gpa_dir) in sys.path:
            sys.path.remove(str(gpa_dir))

    if not result.get("success"):
        return result

    # Normalize summary fields
    results = result.get("results", {})
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
        "report_md": result.get("report_md", ""),
        "raw": result,
    }


def save_gpa_json(result: dict, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, ensure_ascii=False, default=str)
    return path
