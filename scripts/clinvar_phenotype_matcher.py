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
import shutil
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


def _detect_clinvar_chrom_prefix(vcf_path: Path) -> str:
    """Return 'chr' if ClinVar VCF contigs use chr prefix, else ''."""
    result = subprocess.run(
        ["bcftools", "view", "-h", str(vcf_path)],
        capture_output=True, text=True, check=False,
    )
    for line in result.stdout.split("\n"):
        if line.startswith("##contig=<ID="):
            c = line.split("##contig=<ID=")[1].split(",")[0]
            return "chr" if c.startswith("chr") else ""
    return ""


def _normalize_variant_key(chrom: str, pos: int, ref: str, alt: str) -> str:
    """Canonical key for variant lookup (chr stripped, upper case)."""
    return f"{str(chrom).lstrip('chr')}:{pos}:{str(ref).upper()}:{str(alt).upper()}"


def _extract_clnsig(info: str) -> list[str]:
    """Parse all CLNSIG values from ClinVar INFO (pipe-separated)."""
    sigs: list[str] = []
    m = re.search(r"CLNSIG=([^;]+)", info)
    if m:
        for part in m.group(1).split("|"):
            part = part.strip()
            if part and part not in sigs:
                sigs.append(part)
    return sigs


def _extract_review_status(info: str) -> str:
    """Parse CLNREVSTAT (e.g. 'criteria_provided', 'single_submitter')."""
    m = re.search(r"CLNREVSTAT=([^;]+)", info)
    if m:
        return m.group(1).replace("_", " ").strip()
    return ""


def _extract_clinvar_accession(info: str) -> str:
    """Parse first ClinVar Variation ID from ID field or INFO."""
    m = re.search(r"ALLELEID=([^;]+)", info)
    if m:
        return m.group(1).split("|")[0].strip()
    return ""


def _batch_annotate_clinvar(
    variants: list[dict],
    vcf_path: Path = CLINVAR_VCF,
) -> dict[str, dict]:
    """Batch-annotate a list of variants from the local ClinVar VCF.

    Returns a dict mapping normalized variant keys to annotation dicts with
    clinvar_diseases, clinvar_sig, clinvar_sigs (all), clinvar_review_status,
    and clinvar_accession.
    """
    annotations: dict[str, dict] = {}
    if not vcf_path.exists() or not variants:
        return annotations

    prefix = _detect_clinvar_chrom_prefix(vcf_path)

    # Build a BED of positions to query (0-based half-open)
    tmp_dir = Path(tempfile.mkdtemp(prefix="drq_clinvar_batch_"))
    try:
        bed_path = tmp_dir / "query.bed"
        with open(bed_path, "w", encoding="utf-8") as f:
            for v in variants:
                chrom = v.get("chrom") or v.get("CHROM")
                pos = v.get("pos") or v.get("POS")
                if not chrom or not pos:
                    continue
                c = f"{prefix}{str(chrom).lstrip('chr')}"
                p = int(pos)
                f.write(f"{c}\t{p - 1}\t{p}\n")

        proc = subprocess.run(
            ["bcftools", "view", "-H", "-R", str(bed_path), str(vcf_path)],
            capture_output=True, text=True, check=False,
        )
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    for line in proc.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 8:
            continue
        chrom, pos, _id, ref, alt, qual, filt, info = parts[:8]
        alts = alt.split(",")
        diseases = _extract_clndn(info)
        sigs = _extract_clnsig(info)
        sig = sigs[0] if sigs else ""
        review = _extract_review_status(info)
        accession = _extract_clinvar_accession(info)
        for a in alts:
            key = _normalize_variant_key(chrom, int(pos), ref, a)
            annotations[key] = {
                "clinvar_diseases": sorted(set(diseases)),
                "clinvar_sig": sig,
                "clinvar_sigs": sigs,
                "clinvar_review_status": review,
                "clinvar_accession": accession,
            }
    return annotations


def filter_variants_by_clinvar_disease(
    variants: list[dict],
    disease_name: str,
    require_match: bool = False,
) -> list[dict]:
    """Tag variants with ClinVar disease match and optionally drop non-matches.

    Each variant is expected to have chrom/pos/ref/alt and optionally
    clinvar_diseases/clinvar_sig. Missing annotations are fetched from the
    local ClinVar VCF in a single batch query.
    """
    keywords = disease_keywords(disease_name)

    # Identify variants that need ClinVar annotation
    to_annotate = [
        v for v in variants
        if v.get("clinvar_diseases") is None
        and (v.get("chrom") or v.get("CHROM"))
        and (v.get("pos") or v.get("POS"))
        and (v.get("ref") or v.get("REF"))
        and (v.get("alt") or v.get("ALT"))
    ]
    annotation_lookup = _batch_annotate_clinvar(to_annotate)

    out = []
    for v in variants:
        chrom = v.get("chrom") or v.get("CHROM")
        pos = v.get("pos") or v.get("POS")
        ref = v.get("ref") or v.get("REF")
        alt = v.get("alt") or v.get("ALT")

        diseases = v.get("clinvar_diseases")
        if diseases is None and chrom and pos and ref and alt:
            key = _normalize_variant_key(chrom, int(pos), ref, alt)
            ann = annotation_lookup.get(key, {"clinvar_diseases": [], "clinvar_sig": ""})
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


def classify_clinvar_sig(sig: str) -> str:
    """Classify a ClinVar significance string into a canonical category.

    Categories: pathogenic, likely_pathogenic, vus, benign, likely_benign,
    risk_factor, protective, drug_response, conflicting, other.
    """
    if not sig:
        return "other"
    s = sig.lower()
    # Conflicting takes precedence
    if "conflicting" in s:
        return "conflicting"
    # Pathogenic / likely pathogenic
    if "pathogenic" in s and "likely_pathogenic" not in s:
        return "pathogenic"
    if "likely_pathogenic" in s:
        return "likely_pathogenic"
    # Benign
    if "likely_benign" in s:
        return "likely_benign"
    if "benign" in s:
        return "benign"
    # VUS
    if "uncertain" in s or "vus" in s or "significance" in s:
        return "vus"
    # Others
    if "risk_factor" in s:
        return "risk_factor"
    if "protective" in s:
        return "protective"
    if "drug_response" in s:
        return "drug_response"
    return "other"


def is_clinvar_plp(sig: str) -> bool:
    """Return True if ClinVar sig is pathogenic or likely pathogenic."""
    cat = classify_clinvar_sig(sig)
    return cat in ("pathogenic", "likely_pathogenic")


def is_clinvar_benign(sig: str) -> bool:
    """Return True if ClinVar sig is benign or likely benign."""
    cat = classify_clinvar_sig(sig)
    return cat in ("benign", "likely_benign")


def is_clinvar_vus(sig: str) -> bool:
    """Return True if ClinVar sig is VUS."""
    return classify_clinvar_sig(sig) == "vus"


def enrich_variants_with_clinvar(
    variants: list[dict],
    disease_name: str,
    vcf_path: Path = CLINVAR_VCF,
) -> list[dict]:
    """Annotate variants with ClinVar information without filtering or downweighting.

    This is the preferred enrichment path for disease-space variants: every
    variant receives its ClinVar classification and disease-match flag, but no
    variant is dropped because of VUS, conflicting, or non-pathogenic status.

    Returns a new list of variants with added keys:
      - clinvar_diseases
      - clinvar_sig, clinvar_sigs
      - clinvar_review_status
      - clinvar_accession
      - clinvar_category (canonical category)
      - clinvar_plp (bool)
      - clinvar_disease_match (bool, against query disease)
    """
    keywords = disease_keywords(disease_name)

    to_annotate = [
        v for v in variants
        if v.get("clinvar_diseases") is None
        and (v.get("chrom") or v.get("CHROM"))
        and (v.get("pos") or v.get("POS"))
        and (v.get("ref") or v.get("REF"))
        and (v.get("alt") or v.get("ALT"))
    ]
    annotation_lookup = _batch_annotate_clinvar(to_annotate, vcf_path=vcf_path)

    out = []
    for v in variants:
        chrom = v.get("chrom") or v.get("CHROM")
        pos = v.get("pos") or v.get("POS")
        ref = v.get("ref") or v.get("REF")
        alt = v.get("alt") or v.get("ALT")

        diseases = v.get("clinvar_diseases")
        ann: dict = {}
        if diseases is None and chrom and pos and ref and alt:
            key = _normalize_variant_key(chrom, int(pos), ref, alt)
            ann = annotation_lookup.get(
                key,
                {
                    "clinvar_diseases": [],
                    "clinvar_sig": "",
                    "clinvar_sigs": [],
                    "clinvar_review_status": "",
                    "clinvar_accession": "",
                },
            )
            v = {**v, **ann}
            diseases = ann["clinvar_diseases"]

        sig = v.get("clinvar_sig", "")
        v["clinvar_category"] = classify_clinvar_sig(sig)
        v["clinvar_plp"] = is_clinvar_plp(sig)

        # Disease keyword match (informational only)
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
    # Pre-filter to P/LP records only; this avoids parsing the entire ClinVar VCF.
    proc = subprocess.run(
        [
            "bcftools", "view", "-H",
            '-i', 'CLNSIG~"pathogenic" || CLNSIG~"likely_pathogenic"',
            str(vcf_path),
        ],
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
            if not _clinvar_matches_disease(info, keywords):
                continue
            start = int(pos) - 1
            end = start + len(ref)
            sig = "pathogenic"
            m = re.search(r"CLNSIG=([^;]+)", info)
            if m:
                sig = m.group(1).lower()
            for a in alt.split(","):
                f.write(f"{chrom}\t{start}\t{end}\t{ref}>{a}\t{sig}\n")
                count += 1
    logger.info("Wrote %d disease-matched ClinVar variants to %s", count, output_bed)
    return output_bed
