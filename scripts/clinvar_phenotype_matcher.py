"""Match ClinVar variants to disease phenotypes using the local ClinVar VCF.

The dgra-prefilter ClinVar BED used as a safety net only contains coordinates.
To avoid pulling in pathogenic variants for unrelated diseases, we annotate each
variant with ClinVar's CLNDN/CLNDNINCL disease names and keep only those that
match the query disease keywords.
"""

from __future__ import annotations

import gzip
import logging
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Iterable, Optional

from constants import CLINVAR_VCF, DISEASE_CLINVAR_KEYWORDS, resolve_builtin_disease_key

logger = logging.getLogger(__name__)


def _normalize_text(text: str) -> str:
    """Lowercase, remove punctuation, collapse whitespace."""
    return re.sub(r"[^a-z0-9 ]+", " ", text.lower()).strip()


def _extract_clndn(info: str) -> list[str]:
    """Parse CLNDN and CLNDNINCL fields from ClinVar INFO."""
    names: list[str] = []
    for field in ("CLNDN", "CLNDNINCL"):
        # Match either CLNDN=... or CLNDN=...;
        m = re.search(rf"{field}=([^;]+)", info)
        if m:
            val = m.group(1)
            # ClinVar uses | and , as separators, and "not_provided"
            for part in re.split(r"[|,]", val):
                part = part.strip()
                if part and part.lower() != "not_provided":
                    names.append(part)
    return names


def _clinvar_matches_disease(info: str, keywords: Iterable[str]) -> bool:
    """Return True if any ClinVar disease name matches a query keyword."""
    names = _extract_clndn(info)
    if not names:
        return False
    norm_names = [_normalize_text(n) for n in names]
    for kw in keywords:
        kw_norm = _normalize_text(kw)
        if not kw_norm:
            continue
        for name in norm_names:
            if kw_norm in name or name in kw_norm:
                return True
    return False


def _read_clinvar_records(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    vcf_path: Path = CLINVAR_VCF,
) -> list[dict]:
    """Return ClinVar records matching a variant position/alleles."""
    if not vcf_path.exists():
        return []

    # Ensure consistent chromosome style
    result = subprocess.run(
        ["bcftools", "view", "-h", str(vcf_path)],
        capture_output=True, text=True, check=False,
    )
    prefix = ""
    for line in result.stdout.split("\n"):
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
    matches = []
    for line in proc.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 8:
            continue
        if parts[3] != ref:
            continue
        alts = parts[4].split(",")
        if alt not in alts:
            continue
        matches.append({
            "chrom": parts[0],
            "pos": int(parts[1]),
            "ref": parts[3],
            "alt": alt,
            "info": parts[7],
        })
    return matches


def annotate_variant_clinvar_disease(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    vcf_path: Path = CLINVAR_VCF,
) -> dict:
    """Annotate a variant with ClinVar disease names and match status.

    Returns dict with:
      - clinvar_diseases: list of disease names
      - clinvar_sig: parsed CLNSIG (first value)
      - match_disease: bool (whether any disease name matches query keywords)
    """
    records = _read_clinvar_records(chrom, pos, ref, alt, vcf_path)
    diseases: list[str] = []
    sigs: list[str] = []
    for rec in records:
        diseases.extend(_extract_clndn(rec["info"]))
        m = re.search(r"CLNSIG=([^;]+)", rec["info"])
        if m:
            for part in m.group(1).split("|"):
                part = part.strip()
                if part and part not in sigs:
                    sigs.append(part)

    diseases = sorted(set(diseases))
    sig = sigs[0] if sigs else ""
    return {
        "clinvar_diseases": diseases,
        "clinvar_sig": sig,
    }


def disease_keywords(disease_name: str) -> list[str]:
    """Return ClinVar phenotype keywords for a disease query."""
    # Use canonical English key for keyword lookup so Chinese queries map to
    # the correct keyword set.
    canonical = resolve_builtin_disease_key(disease_name)
    lookup_name = canonical or disease_name
    norm = _normalize_text(lookup_name)
    keywords = set(DISEASE_CLINVAR_KEYWORDS.get(norm, []))
    keywords.add(lookup_name)
    # Add disease stem (e.g. "alzheimer disease" -> "alzheimer")
    if " " in lookup_name:
        keywords.add(lookup_name.split()[0])
    return sorted(keywords)


def filter_variants_by_clinvar_disease(
    variants: list[dict],
    disease_name: str,
    require_match: bool = False,
) -> list[dict]:
    """Tag variants with ClinVar disease match and optionally drop non-matches.

    Each variant is expected to have chrom/pos/ref/alt and optionally
    clinvar_diseases/clinvar_sig. Missing annotations are fetched from the
    local ClinVar VCF on demand.
    """
    keywords = disease_keywords(disease_name)
    out = []
    for v in variants:
        chrom = v.get("chrom") or v.get("CHROM")
        pos = v.get("pos") or v.get("POS")
        ref = v.get("ref") or v.get("REF")
        alt = v.get("alt") or v.get("ALT")

        diseases = v.get("clinvar_diseases")
        if diseases is None and chrom and pos and ref and alt:
            ann = annotate_variant_clinvar_disease(str(chrom), int(pos), str(ref), str(alt))
            v = {**v, **ann}
            diseases = ann["clinvar_diseases"]

        norm_diseases = [_normalize_text(d) for d in (diseases or [])]
        matched = False
        for kw in keywords:
            kw_norm = _normalize_text(kw)
            for name in norm_diseases:
                if kw_norm in name or name in kw_norm:
                    matched = True
                    break
            if matched:
                break

        v["clinvar_disease_match"] = matched
        if not require_match or matched:
            out.append(v)
    return out


def build_disease_clinvar_bed(
    disease_name: str,
    output_bed: Path,
    vcf_path: Path = CLINVAR_VCF,
) -> Path:
    """Create a BED of ClinVar P/LP variants whose phenotype matches disease.

    Used to replace the broad ClinVar pathogenic safety net with a disease-
    relevant one.
    """
    keywords = disease_keywords(disease_name)
    output_bed.parent.mkdir(parents=True, exist_ok=True)

    logger.info(
        "Building disease-matched ClinVar BED for '%s' from %s",
        disease_name, vcf_path,
    )
    proc = subprocess.run(
        ["bcftools", "view", "-H", str(vcf_path)],
        capture_output=True, text=True, check=False,
    )
    count = 0
    with open(output_bed, "w", encoding="utf-8") as f:
        for line in proc.stdout.strip().split("\n"):
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _id, ref, alt, qual, filt, info = parts[:8]
            # Keep pathogenic/likely_pathogenic only
            m = re.search(r"CLNSIG=([^;]+)", info)
            if not m:
                continue
            sig = m.group(1).lower()
            if not any(s in sig for s in ("pathogenic", "likely_pathogenic")):
                continue
            if not _clinvar_matches_disease(info, keywords):
                continue
            start = int(pos) - 1
            end = start + len(ref)
            for a in alt.split(","):
                f.write(f"{chrom}\t{start}\t{end}\t{ref}>{a}\t{sig}\n")
                count += 1
    logger.info("Wrote %d disease-matched ClinVar variants to %s", count, output_bed)
    return output_bed
