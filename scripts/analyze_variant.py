#!/usr/bin/env python3
"""
GRCh38 Variant Functional Impact Analyzer

Analyzes a single GRCh38 variant by querying public APIs and local databases:
- Ensembl VEP (SIFT, PolyPhen, AlphaMissense, CADD, SpliceAI, frequencies)
- gnomAD v4.1 GraphQL (population AF)
- NCBI ClinVar (classification, phenotypes)
- Local OMIM SQLite (gene-phenotype associations)
- UniProt (protein function and features)
- Europe PMC (variant-specific literature)

Output: structured JSON and/or Markdown summary.
"""

from __future__ import annotations

import argparse
import json
import re
import sqlite3
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional
from urllib.parse import urlencode

import requests

DEFAULT_OMIM_DB = Path.home() / ".workbuddy" / "data" / "omim" / "omim.db"


class Variant:
    """Canonical GRCh38 variant representation."""

    def __init__(
        self,
        raw: str,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        rsid: Optional[str] = None,
        hgvs_g: Optional[str] = None,
    ):
        self.raw = raw
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.rsid = rsid
        self.hgvs_g = hgvs_g or f"{chrom}:g.{pos}{ref}>{alt}"

    def __repr__(self) -> str:
        return f"Variant({self.hgvs_g})"


def _normalize_chrom(chrom: str) -> str:
    chrom = chrom.strip().lower().replace("chr", "")
    if chrom in ("x", "y", "m", "mt"):
        return chrom.upper()
    return chrom


def parse_variant(variant_str: str) -> Variant:
    """Parse a variant from multiple supported formats."""

    variant_str = variant_str.strip()

    # rsID
    if variant_str.lower().startswith("rs"):
        return Variant(
            raw=variant_str,
            chrom="",
            pos=0,
            ref="",
            alt="",
            rsid=variant_str.lower(),
            hgvs_g="",
        )

    # HGVS genomic: 2:g.21007456G>C
    m = re.match(r"^(chr)?([0-9XYMTxymt]+):g\.(\d+)([ACGTNacgtn]+)>([ACGTNacgtn]+)$", variant_str)
    if m:
        chrom = _normalize_chrom(m.group(2))
        pos = int(m.group(3))
        ref = m.group(4).upper()
        alt = m.group(5).upper()
        return Variant(
            raw=variant_str,
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            hgvs_g=f"{chrom}:g.{pos}{ref}>{alt}",
        )

    # chrom:pos:ref:alt
    parts = variant_str.replace(":", " ").split()
    if len(parts) == 4:
        chrom = _normalize_chrom(parts[0])
        pos = int(parts[1])
        ref = parts[2].upper()
        alt = parts[3].upper()
        return Variant(
            raw=variant_str,
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            hgvs_g=f"{chrom}:g.{pos}{ref}>{alt}",
        )

    # HGVS coding (e.g., NM_000384.3:c.9412C>G): leave resolution to VEP
    m = re.match(r"^(NM_[0-9]+\.[0-9]+):c\.(.+)$", variant_str)
    if m:
        return Variant(
            raw=variant_str,
            chrom="",
            pos=0,
            ref="",
            alt="",
            hgvs_g=variant_str,
        )

    raise ValueError(
        f"Unsupported variant format: {variant_str}. "
        "Expected one of: chr2:21007456:G:C, 2:g.21007456G>C, rs755753065, "
        "NM_000384.3:c.9412C>G"
    )


def _session(proxy: Optional[str] = None) -> requests.Session:
    """Create a requests session with explicit proxy handling."""
    sess = requests.Session()
    sess.trust_env = False
    if proxy:
        sess.proxies = {"http": proxy, "https": proxy}
    return sess


def _request(
    sess: requests.Session,
    method: str,
    url: str,
    **kwargs,
) -> Optional[Any]:
    """Make an HTTP request with retries and return JSON."""
    for attempt in range(3):
        try:
            r = sess.request(method, url, timeout=60, **kwargs)
            if r.status_code == 200:
                return r.json()
            if r.status_code in (429, 502, 503, 504):
                time.sleep(2 ** attempt)
                continue
            return None
        except Exception:
            if attempt == 2:
                return None
            time.sleep(2 ** attempt)
    return None


# ---------------------------------------------------------------------------
# Ensembl VEP
# ---------------------------------------------------------------------------


def resolve_rsid(rsid: str, proxy: Optional[str] = None) -> Optional[Variant]:
    """Resolve an rsID to chrom/pos/ref/alt using Ensembl variation API."""
    sess = _session(proxy)
    url = f"https://rest.ensembl.org/variation/human/{rsid.lower()}"
    data = _request(sess, "GET", url, headers={"Content-Type": "application/json"}, params={"pops": 0})
    if not data or "mappings" not in data:
        return None
    for mapping in data["mappings"]:
        if mapping.get("assembly_name") == "GRCh38":
            chrom = _normalize_chrom(mapping.get("seq_region_name", ""))
            pos = mapping.get("start")
            allele_string = mapping.get("allele_string", "")
            parts = allele_string.split("/") if allele_string else []
            if chrom and pos and len(parts) == 2:
                return Variant(
                    raw=rsid,
                    chrom=chrom,
                    pos=int(pos),
                    ref=parts[0].upper(),
                    alt=parts[1].upper(),
                    rsid=rsid.lower(),
                    hgvs_g=f"{chrom}:g.{pos}{parts[0].upper()}>{parts[1].upper()}",
                )
    return None


def query_vep(variant: Variant, proxy: Optional[str] = None) -> tuple[Dict[str, Any], Variant]:
    """Query Ensembl VEP with AlphaMissense, CADD, and SpliceAI.

    Returns a tuple of (vep_result_dict, updated_variant).
    """
    sess = _session(proxy)
    url = "https://rest.ensembl.org/vep/human/hgvs"
    payload = {
        "hgvs_notations": [variant.hgvs_g],
        "AlphaMissense": True,
        "CADD": True,
        "SpliceAI": True,
    }
    data = _request(sess, "POST", url, headers={"Content-Type": "application/json"}, json=payload)
    if not data or not isinstance(data, list):
        return {"error": "VEP returned no data"}, variant

    top = data[0]
    updated_variant = _update_variant_from_vep(variant, top)
    transcript = _pick_canonical_transcript(top.get("transcript_consequences", []))

    rsid = _resolve_rsid_from_vep(top, updated_variant.rsid)
    cdna, protein = _extract_hgvs(transcript)

    result: Dict[str, Any] = {
        "query": updated_variant.hgvs_g,
        "input": updated_variant.raw,
        "rsid": rsid,
        "gene_symbol": transcript.get("gene_symbol") if transcript else None,
        "transcript": transcript.get("transcript_id") if transcript else None,
        "cdna": cdna,
        "protein": protein,
        "protein_start": transcript.get("protein_start") if transcript else None,
        "protein_end": transcript.get("protein_end") if transcript else None,
        "amino_acids": transcript.get("amino_acids") if transcript else None,
        "consequence_terms": transcript.get("consequence_terms") if transcript else None,
        "sift": _extract_sift(transcript),
        "polyphen": _extract_polyphen(transcript),
        "alphamissense": _extract_alphamissense(transcript),
        "cadd_phred": _extract_cadd(top, transcript),
        "spliceai": _extract_spliceai(transcript),
        "gnomad_frequencies": _extract_frequencies(top),
        "all_transcript_consequences": top.get("transcript_consequences", []),
    }
    return result, updated_variant


def _pick_canonical_transcript(consequences: List[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    if not consequences:
        return None
    # Prefer MANE select, then canonical, then first with most severe consequence
    for tc in consequences:
        if tc.get("mane_select"):
            return tc
    for tc in consequences:
        if tc.get("canonical") == 1:
            return tc
    return consequences[0]


def _resolve_rsid_from_vep(top: Dict[str, Any], fallback: Optional[str]) -> Optional[str]:
    if fallback and fallback.lower().startswith("rs"):
        return fallback
    for cv in top.get("colocated_variants", []):
        cid = cv.get("id", "")
        if isinstance(cid, str) and cid.lower().startswith("rs"):
            return cid
    return fallback


def _update_variant_from_vep(variant: Variant, top: Dict[str, Any]) -> Variant:
    """Resolve chrom/pos/ref/alt/rsid from VEP response for rsID/HGVS coding inputs."""
    rsid = _resolve_rsid_from_vep(top, variant.rsid)

    if variant.chrom and variant.pos and variant.ref and variant.alt:
        if rsid != variant.rsid:
            return Variant(
                raw=variant.raw,
                chrom=variant.chrom,
                pos=variant.pos,
                ref=variant.ref,
                alt=variant.alt,
                rsid=rsid,
                hgvs_g=variant.hgvs_g,
            )
        return variant

    chrom = top.get("seq_region_name")
    pos = top.get("start")
    allele_string = top.get("allele_string", "")
    parts = allele_string.split("/") if allele_string else []
    if chrom and pos and len(parts) == 2:
        return Variant(
            raw=variant.raw,
            chrom=_normalize_chrom(str(chrom)),
            pos=int(pos),
            ref=parts[0].upper(),
            alt=parts[1].upper(),
            rsid=rsid,
            hgvs_g=f"{_normalize_chrom(str(chrom))}:g.{pos}{parts[0].upper()}>{parts[1].upper()}",
        )
    return variant


def _extract_sift(tc: Optional[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    if not tc:
        return None
    pred = tc.get("sift_prediction")
    score = tc.get("sift_score")
    if pred is None and score is None:
        return None
    return {"prediction": pred, "score": score}


def _extract_polyphen(tc: Optional[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    if not tc:
        return None
    pred = tc.get("polyphen_prediction")
    score = tc.get("polyphen_score")
    if pred is None and score is None:
        return None
    return {"prediction": pred, "score": score}


def _extract_alphamissense(tc: Optional[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    if not tc:
        return None
    am = tc.get("alphamissense", {})
    if not am:
        return None
    return {
        "class": am.get("am_class"),
        "pathogenicity": am.get("am_pathogenicity"),
    }


def _extract_cadd(top: Dict[str, Any], transcript: Optional[Dict[str, Any]]) -> Optional[float]:
    # VEP returns CADD in transcript consequence; top-level may not have it.
    if transcript and transcript.get("cadd_phred") is not None:
        return transcript.get("cadd_phred")
    return top.get("cadd_phred")


def _short_protein_change(hgvs_p: str) -> Optional[str]:
    """Convert p.Pro3138Ala -> p.P3138A for literature search."""
    m = re.match(r"^p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})$", hgvs_p)
    if not m:
        return None
    ref = m.group(1)[0].upper()
    pos = m.group(2)
    alt = m.group(3)[0].upper()
    return f"p.{ref}{pos}{alt}"


def _aa_long(aa: str) -> str:
    """Convert single-letter amino acid to 3-letter abbreviation."""
    table = {
        "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
        "E": "Glu", "Q": "Gln", "G": "Gly", "H": "His", "I": "Ile",
        "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
        "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
        "*": "Ter", "X": "Xaa",
    }
    return table.get(aa.upper(), aa)


def _extract_hgvs(tc: Optional[Dict[str, Any]]) -> tuple[Optional[str], Optional[str]]:
    """Extract or construct cDNA and protein HGVS from VEP transcript consequence."""
    if not tc:
        return None, None

    cdna = tc.get("hgvsc")
    protein = tc.get("hgvsp")
    if cdna and protein:
        return cdna, protein

    # Construct from available fields
    codons = tc.get("codons", "")
    aa = tc.get("amino_acids", "")
    cds_start = tc.get("cds_start")
    protein_start = tc.get("protein_start")

    if not cdna and codons and cds_start:
        ref, alt = codons.split("/") if "/" in codons else (codons, "")
        if len(ref) == 3 and len(alt) == 3:
            ref_base = ref[0].upper()
            alt_base = alt[0].upper()
            cdna = f"c.{cds_start}{ref_base}>{alt_base}"

    if not protein and aa and protein_start:
        ref_aa, alt_aa = aa.split("/") if "/" in aa else (aa, "")
        if len(ref_aa) == 1 and len(alt_aa) == 1:
            protein = f"p.{_aa_long(ref_aa)}{protein_start}{_aa_long(alt_aa)}"

    return cdna, protein


def _extract_spliceai(tc: Optional[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    if not tc:
        return None
    sai = tc.get("spliceai", {})
    if not sai:
        return None
    return {
        "DS_AG": sai.get("DS_AG"),
        "DS_AL": sai.get("DS_AL"),
        "DS_DG": sai.get("DS_DG"),
        "DS_DL": sai.get("DS_DL"),
    }


def _extract_frequencies(top: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    colocated = top.get("colocated_variants", [])
    if not colocated:
        return None
    freqs = colocated[0].get("frequencies", {})
    alt_freqs = next(iter(freqs.values()), {}) if freqs else {}
    if not alt_freqs:
        return None
    return {
        k: v
        for k, v in alt_freqs.items()
        if k.startswith("gnomad") and v is not None
    } or None


# ---------------------------------------------------------------------------
# gnomAD
# ---------------------------------------------------------------------------


def query_gnomad(variant: Variant, proxy: Optional[str] = None) -> Dict[str, Any]:
    """Query gnomAD v4.1 GraphQL for exome/genome frequencies."""
    if not variant.chrom or not variant.pos:
        return {"error": "Cannot query gnomAD without chrom/pos"}

    sess = _session(proxy)
    query = """
    query GetVariant($chrom: String!, $start: Int!, $stop: Int!, $dataset: DatasetId!) {
      region(chrom: $chrom, start: $start, stop: $stop, reference_genome: GRCh38) {
        variants(dataset: $dataset) {
          variant_id
          exome { ac an af }
          genome { ac an af }
        }
      }
    }
    """
    variables = {
        "chrom": variant.chrom,
        "start": max(1, variant.pos - 5),
        "stop": variant.pos + 5,
        "dataset": "gnomad_r4",
    }
    payload = {"query": query, "variables": variables}
    data = _request(
        sess,
        "POST",
        "https://gnomad.broadinstitute.org/api/",
        json=payload,
        headers={"Content-Type": "application/json"},
    )
    if not data or "data" not in data:
        return {"error": "gnomAD returned no data"}

    target_id = f"{variant.chrom}-{variant.pos}-{variant.ref}-{variant.alt}"
    variants = data["data"]["region"]["variants"]
    for v in variants:
        if v["variant_id"] == target_id:
            return {
                "variant_id": v["variant_id"],
                "exome": v.get("exome"),
                "genome": v.get("genome"),
            }
    return {"error": f"Variant {target_id} not found in gnomAD"}


# ---------------------------------------------------------------------------
# ClinVar
# ---------------------------------------------------------------------------


def query_clinvar(variant: Variant, vep: Dict[str, Any], proxy: Optional[str] = None) -> Dict[str, Any]:
    """Query ClinVar by rsID or HGVS."""
    sess = _session(proxy)

    # Try rsID via esearch -> esummary
    rsid = variant.rsid or vep.get("rsid")
    if rsid and rsid.lower().startswith("rs"):
        uid = _clinvar_search_uid(sess, rsid)
        if uid:
            data = _clinvar_esummary(sess, uid)
            if data:
                return data

    # Fallback: search by HGVS
    if variant.hgvs_g:
        uid = _clinvar_search_uid(sess, variant.hgvs_g)
        if uid:
            data = _clinvar_esummary(sess, uid)
            if data:
                return data

    return {"error": "No ClinVar annotation found"}


def _clinvar_search_uid(sess: requests.Session, term: str) -> Optional[str]:
    """Return the first ClinVar UID matching a search term."""
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {"db": "clinvar", "term": term, "retmode": "json", "retmax": 1}
    data = _request(sess, "GET", url, params=params)
    if not data or "esearchresult" not in data:
        return None
    ids = data["esearchresult"].get("idlist", [])
    return ids[0] if ids else None


def _clinvar_esummary(sess: requests.Session, uid: str) -> Optional[Dict[str, Any]]:
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    params = {"db": "clinvar", "id": uid, "retmode": "json"}
    data = _request(sess, "GET", url, params=params)
    if not data or "result" not in data:
        return None
    result_obj = data["result"]
    uids = result_obj.get("uids", [])
    if not uids:
        return None
    raw = result_obj[uids[0]]
    return _parse_clinvar_record(raw)



def _parse_clinvar_record(raw: Dict[str, Any]) -> Dict[str, Any]:
    germline = raw.get("germline_classification", {})
    traits = []
    # The trait_set is nested inside germline_classification in current ClinVar esummary.
    for trait_set in germline.get("trait_set") or raw.get("trait_set") or []:
        if isinstance(trait_set, dict):
            name = trait_set.get("trait_name")
            if isinstance(name, str):
                traits.append(name)
            elif isinstance(name, list):
                traits.extend(str(n) for n in name)
    return {
        "accession": raw.get("accession"),
        "accession_version": raw.get("accession_version"),
        "title": raw.get("title"),
        "classification": germline.get("description"),
        "review_status": germline.get("review_status"),
        "last_evaluated": germline.get("last_evaluated"),
        "traits": traits,
    }


# ---------------------------------------------------------------------------
# OMIM
# ---------------------------------------------------------------------------


def query_omim(vep: Dict[str, Any], omim_db: Optional[Path] = None) -> Optional[Dict[str, Any]]:
    """Query local OMIM SQLite for gene-phenotype associations."""
    db_path = omim_db or DEFAULT_OMIM_DB
    if not db_path.exists():
        return None

    gene_symbol = _extract_gene_symbol(vep)
    if not gene_symbol:
        return None

    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    c = conn.cursor()

    # Try exact HGNC approved symbol (column may contain leading/trailing spaces)
    c.execute(
        "SELECT mim_number, title, geneMap FROM omim WHERE TRIM(hgnc_approved_gene_symbol) = ?",
        (gene_symbol,),
    )
    row = c.fetchone()
    if not row:
        # Search title/symbols, but avoid substring matches (e.g., APOB in HNRNPAB)
        c.execute(
            "SELECT mim_number, title, geneMap FROM omim WHERE "
            "title LIKE ? OR title LIKE ? OR TRIM(symbols) = ? OR symbols LIKE ? OR symbols LIKE ?",
            (
                f"%{gene_symbol};%",
                f"%{gene_symbol} (%",
                gene_symbol,
                f"%{gene_symbol}\n%",
                f"{gene_symbol};%",
            ),
        )
        row = c.fetchone()
    conn.close()

    if not row:
        return None

    gene_map = _safe_json(row["geneMap"])
    phenotypes = _parse_omim_phenotypes(gene_map)

    return {
        "mim_number": row["mim_number"],
        "title": row["title"],
        "gene_symbol": gene_symbol,
        "phenotypes": phenotypes,
    }


def _extract_gene_symbol(vep: Dict[str, Any]) -> Optional[str]:
    return vep.get("gene_symbol")


def _safe_json(value: Any) -> Any:
    if value is None:
        return None
    if isinstance(value, str):
        try:
            return json.loads(value)
        except json.JSONDecodeError:
            return value
    return value


def _parse_omim_phenotypes(gene_map: Any) -> List[Dict[str, Any]]:
    """Parse OMIM geneMap/phenotypeMap into a normalized list."""
    phenotypes: List[Dict[str, Any]] = []
    if not isinstance(gene_map, list):
        return phenotypes
    for entry in gene_map:
        if not isinstance(entry, dict):
            continue
        # geneMap format
        name = entry.get("Phenotype View")
        if name:
            phenotypes.append({
                "phenotype": name,
                "mim_number": entry.get("Phenotype MIM number"),
                "inheritance": entry.get("Inheritance"),
            })
            continue
        # phenotypeMap / nested format
        phenos = entry.get("phenotypes")
        if isinstance(phenos, list):
            for p in phenos:
                if isinstance(p, dict):
                    phenotypes.append({
                        "phenotype": p.get("phenotype"),
                        "mim_number": p.get("mim_number"),
                        "inheritance": p.get("inheritance"),
                    })
    return phenotypes


# ---------------------------------------------------------------------------
# UniProt
# ---------------------------------------------------------------------------


def query_uniprot(gene_symbol: str, protein_pos: Optional[int] = None, proxy: Optional[str] = None) -> Dict[str, Any]:
    """Query UniProt for protein function and features near the variant."""
    sess = _session(proxy)

    # First map gene symbol to UniProt accession
    search_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"gene:{gene_symbol} AND organism_id:9606 AND reviewed:true",
        "fields": "accession,id,gene_names,protein_name,length",
        "format": "json",
        "size": 1,
    }
    search_data = _request(sess, "GET", search_url, params=params)
    if not search_data or not search_data.get("results"):
        return {"error": f"No UniProt entry found for {gene_symbol}"}

    entry = search_data["results"][0]
    accession = entry.get("primaryAccession")

    # Fetch full entry
    full_data = _request(sess, "GET", f"https://rest.uniprot.org/uniprotkb/{accession}.json")
    if not full_data:
        return {"error": f"Could not fetch UniProt entry {accession}"}

    function_texts = []
    for comment in full_data.get("comments", []):
        if comment.get("commentType") == "FUNCTION":
            for text in comment.get("texts", []):
                function_texts.append(text.get("value", ""))

    features_near = []
    if protein_pos:
        for feature in full_data.get("features", []):
            loc = feature.get("location", {})
            start = loc.get("start", {}).get("value")
            end = loc.get("end", {}).get("value")
            if start and end and int(start) <= protein_pos + 2 and int(end) >= protein_pos - 2:
                features_near.append({
                    "type": feature.get("type"),
                    "description": feature.get("description"),
                    "start": start,
                    "end": end,
                })

    return {
        "accession": accession,
        "gene_symbol": gene_symbol,
        "protein_length": full_data.get("sequence", {}).get("length"),
        "protein_name": full_data.get("proteinDescription", {})
        .get("recommendedName", {})
        .get("fullName", {})
        .get("value"),
        "function": " ".join(function_texts),
        "features_near_variant": features_near,
    }


# ---------------------------------------------------------------------------
# Literature
# ---------------------------------------------------------------------------


def query_literature(variant: Variant, vep: Dict[str, Any], proxy: Optional[str] = None) -> Dict[str, Any]:
    """Search Europe PMC for variant-specific literature."""
    sess = _session(proxy)

    # Build query from protein change if available, otherwise HGVS / rsID
    protein = vep.get("protein")
    cdna = vep.get("cdna")
    parts = []
    if protein:
        parts.append(f'"{protein}"')
        short = _short_protein_change(protein)
        if short:
            parts.append(f'"{short}"')
    if cdna and cdna.startswith("c."):
        parts.append(f'"{cdna}"')
    if variant.rsid:
        parts.append(f'"{variant.rsid}"')
    elif variant.hgvs_g and not variant.hgvs_g.startswith("NM_"):
        parts.append(f'"{variant.hgvs_g}"')

    query = " OR ".join(parts) if parts else f'"{variant.raw}"'
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    params = {"query": query, "format": "json", "pageSize": 10}
    data = _request(sess, "GET", url, params=params)
    if not data:
        return {"query": query, "count": 0, "articles": []}

    articles = []
    for rslt in data.get("resultList", {}).get("result", []):
        articles.append({
            "title": rslt.get("title"),
            "authors": rslt.get("authorString"),
            "journal": rslt.get("journalTitle"),
            "year": rslt.get("pubYear"),
            "pmid": rslt.get("pmid"),
            "doi": rslt.get("doi"),
        })
    return {"query": query, "count": data.get("hitCount", 0), "articles": articles}


# ---------------------------------------------------------------------------
# Interpretation
# ---------------------------------------------------------------------------


def build_interpretation(vep: Dict[str, Any], clinvar: Dict[str, Any]) -> Dict[str, Any]:
    """Generate a concise interpretation summary."""
    evidence: List[str] = []
    clinvar_class = clinvar.get("classification")
    if clinvar_class:
        evidence.append(f"ClinVar: {clinvar_class}")

    am = vep.get("alphamissense") or {}
    am_class = am.get("class")
    if am_class:
        evidence.append(f"AlphaMissense: {am_class}")

    cadd = vep.get("cadd_phred")
    if cadd is not None:
        evidence.append(f"CADD phred {cadd}")

    sift = vep.get("sift") or {}
    if sift.get("prediction"):
        evidence.append(f"SIFT: {sift['prediction']} ({sift.get('score')})")

    polyphen = vep.get("polyphen") or {}
    if polyphen.get("prediction"):
        evidence.append(f"PolyPhen: {polyphen['prediction']} ({polyphen.get('score')})")

    # Summary label
    if clinvar_class and "pathogenic" in clinvar_class.lower():
        summary = "Likely pathogenic (ClinVar)"
    elif clinvar_class and "benign" in clinvar_class.lower():
        summary = "Likely benign (ClinVar)"
    elif am_class == "likely_pathogenic" or (cadd is not None and cadd >= 20):
        summary = "VUS with suspicious in-silico evidence"
    elif am_class == "likely_benign":
        summary = "VUS-low / likely benign"
    else:
        summary = "VUS"

    return {"summary": summary, "evidence": evidence}


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def _val(v: Any, default: str = "N/A") -> str:
    """Format a value for Markdown display, returning default for None/empty."""
    if v is None or v == "" or v == []:
        return default
    return str(v)


def _interpretation_cn(summary_en: str) -> str:
    """Translate interpretation summary to Chinese."""
    mapping = {
        "Likely benign (ClinVar)": "可能良性 (ClinVar)",
        "Pathogenic (ClinVar)": "致病 (ClinVar)",
        "Likely pathogenic (ClinVar)": "可能致病 (ClinVar)",
        "VUS-low / likely benign": "意义未明-低 / 可能良性",
        "VUS with suspicious in-silico evidence": "意义未明 + 计算预测偏高风险",
        "VUS": "意义未明 (VUS)",
    }
    return mapping.get(summary_en, summary_en)


def build_markdown(result: Dict[str, Any]) -> str:
    v = result["variant"]
    vep = result.get("vep", {})
    gnomad = result.get("gnomad", {})
    clinvar = result.get("clinvar", {})
    omim = result.get("omim", {})
    uniprot = result.get("uniprot", {})
    lit = result.get("literature", {})
    interp = result.get("interpretation", {})

    lines = [
        f"# 变异功能影响分析报告：{v['hgvs_g']}",
        "",
        f"**输入**：{v['raw']}  ",
        f"**基因组坐标**：{v['hgvs_g']}  ",
        f"**转录本 / 蛋白**：{_val(vep.get('transcript'))} / {_val(vep.get('protein'))}  ",
        f"**解读**：{_interpretation_cn(interp.get('summary', ''))}",
        "",
        "## 1. 功能预测 (VEP)",
        "",
        "| 工具 | 预测结果 | 分值 |",
        "|------|----------|------|",
    ]

    sift = vep.get("sift") or {}
    polyphen = vep.get("polyphen") or {}
    am = vep.get("alphamissense") or {}
    lines.append(f"| SIFT | {_val(sift.get('prediction'))} | {_val(sift.get('score'))} |")
    lines.append(f"| PolyPhen | {_val(polyphen.get('prediction'))} | {_val(polyphen.get('score'))} |")
    lines.append(f"| AlphaMissense | {_val(am.get('class'))} | {_val(am.get('pathogenicity'))} |")
    lines.append(f"| CADD phred | {_val(vep.get('cadd_phred'))} | — |")

    sai = vep.get("spliceai") or {}
    ds = [sai.get(k) for k in ["DS_AG", "DS_AL", "DS_DG", "DS_DL"]]
    no_effect = all(v == 0 for v in ds if v is not None)
    ds_str = "/".join(str(v if v is not None else "—") for v in ds)
    lines.append(f"| SpliceAI | {'无剪接影响' if no_effect else '可能存在剪接影响'} | DS_AG/AL/DG/DL = {ds_str} |")

    lines.extend(["", "## 2. 人群频率", "", "| 来源 | AC | AN | AF |", "|------|----|----|----|"])
    source_names = {"exome": "外显子组", "genome": "全基因组"}
    for source in ["exome", "genome"]:
        g = gnomad.get(source)
        if isinstance(g, dict):
            lines.append(f"| gnomAD {source_names.get(source, source)} | {_val(g.get('ac'))} | {_val(g.get('an'))} | {_val(g.get('af'))} |")
        else:
            lines.append(f"| gnomAD {source_names.get(source, source)} | — | — | {_val(g)} |")

    vep_freqs = vep.get("gnomad_frequencies") or {}
    freq_display = {
        "gnomade": "VEP gnomAD 外显子组",
        "gnomadg": "VEP gnomAD 全基因组",
        "gnomade_eas": "VEP 东亚外显子组",
        "gnomadg_eas": "VEP 东亚全基因组",
    }
    for k, label in freq_display.items():
        val = vep_freqs.get(k)
        if val is not None:
            lines.append(f"| {label} | — | — | {val} |")
    for k, val in sorted(vep_freqs.items()):
        if k not in freq_display and val:
            lines.append(f"| VEP {k} | — | — | {val} |")

    lines.extend(["", "## 3. ClinVar 注释", ""])
    if clinvar.get("error"):
        lines.append(f"- 未命中：{clinvar['error']}")
    else:
        lines.append(f"- **Accession**：{_val(clinvar.get('accession'))}")
        cls_map = {
            "Pathogenic": "致病",
            "Likely pathogenic": "可能致病",
            "Uncertain significance": "意义未明",
            "Likely benign": "可能良性",
            "Benign": "良性",
        }
        raw_cls = clinvar.get("classification", "")
        chinese_cls = cls_map.get(raw_cls, raw_cls)
        lines.append(f"- **临床意义分类**：{chinese_cls or 'N/A'}（{raw_cls}）")
        review_map = {
            "criteria provided, multiple submitters, no conflicts": "有评估标准，多方提交无冲突",
            "criteria provided, single submitter": "有评估标准，单一提交者",
            "reviewed by expert panel": "专家小组审核",
            "practice guideline": "临床实践指南",
            "no assertion criteria provided": "未提供评估标准",
            "no assertion provided": "无评估声明",
        }
        raw_review = clinvar.get("review_status", "")
        chinese_review = review_map.get(raw_review, raw_review)
        lines.append(f"- **审核状态**：{chinese_review or 'N/A'}")
        last_eval = clinvar.get("last_evaluated")
        if isinstance(last_eval, str):
            last_eval = last_eval.split()[0]
        lines.append(f"- **最后审核日期**：{last_eval or 'N/A'}")
        traits = clinvar.get("traits", [])
        if traits:
            lines.append(f"- **关联表型**：{', '.join(traits)}")
        else:
            lines.append("- **关联表型**：无")

    lines.extend(["", "## 4. OMIM 基因-表型关联", ""])
    if not omim:
        lines.append("- 无本地 OMIM 数据可用。")
    else:
        lines.append(f"- **基因**：{_val(omim.get('gene_symbol'))}（MIM {_val(omim.get('mim_number'))}）")
        lines.append(f"- **条目名称**：{_val(omim.get('title'))}")
        phenos = omim.get("phenotypes", [])
        if phenos:
            lines.append("- **关联疾病 / 表型**：")
            for p in phenos:
                inh_map = {"AD": "常染色体显性", "AR": "常染色体隐性", "XD": "X连锁显性", "XR": "X连锁隐性"}
                inh = inh_map.get(p.get("inheritance", ""), p.get("inheritance", ""))
                lines.append(f"  - {p.get('phenotype')} | MIM {p.get('mim_number')} | {inh}")
        else:
            lines.append("- **关联疾病 / 表型**：无记录")

    lines.extend(["", "## 5. UniProt 蛋白信息", ""])
    if uniprot.get("error"):
        lines.append(f"- 查询失败：{uniprot['error']}")
    else:
        lines.append(f"- **UniProt Accession**：{_val(uniprot.get('accession'))}")
        lines.append(f"- **蛋白名称**：{_val(uniprot.get('protein_name'))}")
        lines.append(f"- **蛋白长度**：{_val(uniprot.get('protein_length'))} 个氨基酸")
        func = uniprot.get("function")
        if func:
            lines.append(f"- **功能描述**：{func}")
        else:
            lines.append("- **功能描述**：无")
        feats = uniprot.get("features_near_variant", [])
        if feats:
            lines.append("- **变异附近的结构特征**：")
            feat_type_map = {
                "Chain": "成熟链",
                "Sequence conflict": "序列冲突",
                "Domain": "结构域",
                "Region": "区域",
                "Binding site": "结合位点",
                "Active site": "活性位点",
                "Modified residue": "修饰残基",
                "Disulfide bond": "二硫键",
            }
            for f in feats:
                chinese_type = feat_type_map.get(f["type"], f["type"])
                lines.append(f"  - {chinese_type}：{f['description']} [{f['start']}-{f['end']}]")
        else:
            lines.append("- **变异附近的结构特征**：无注释")

    lines.extend(["", "## 6. 文献检索", ""])
    lines.append(f"- 检索式：`{lit.get('query')}`")
    lines.append(f"- 命中：{lit.get('count', 0)} 条")
    articles = lit.get("articles", [])
    if articles:
        for a in articles[:5]:
            lines.append(f"  - {a.get('title')}（{a.get('year')}）PMID:{a.get('pmid')}")
    else:
        lines.append("  - 未找到该位点特异性文献。")

    lines.extend(["", "## 7. 综合解读", ""])
    summary_en = interp.get("summary", "")
    summary_cn = _interpretation_cn(summary_en)
    lines.append(f"**{summary_cn}**（{summary_en}）")
    lines.append("")
    lines.append("**证据汇总**：")
    for e in interp.get("evidence", []):
        lines.append(f"- {e}")

    lines.extend(["", "---", "", "*由 grch38-variant-impact skill 生成。本报告仅供研究参考，不可作为临床诊断依据。*"])
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(description="GRCh38 variant functional impact analyzer")
    parser.add_argument("--variant", required=True, help="Variant (chr:pos:ref:alt, HGVS, or rsID)")
    parser.add_argument("--omim-db", type=Path, default=DEFAULT_OMIM_DB, help="Path to OMIM SQLite DB")
    parser.add_argument("--output-dir", type=Path, default=Path("."), help="Output directory")
    parser.add_argument("--format", choices=["json", "md", "both"], default="both", help="Output format")
    parser.add_argument("--proxy", default=None, help="HTTPS proxy URL (e.g., http://127.0.0.1:7890)")
    parser.add_argument("--no-gnomad", action="store_true", help="Skip gnomAD query")
    parser.add_argument("--no-clinvar", action="store_true", help="Skip ClinVar query")
    parser.add_argument("--no-literature", action="store_true", help="Skip literature search")
    parser.add_argument("--no-omim", action="store_true", help="Skip OMIM query")
    args = parser.parse_args()

    try:
        variant = parse_variant(args.variant)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    print(f"Analyzing {variant.raw} ...")

    # Step 0: resolve rsID to coordinates if needed
    if variant.rsid and not variant.chrom:
        resolved = resolve_rsid(variant.rsid, proxy=args.proxy)
        if not resolved:
            print(f"Could not resolve {variant.rsid} to GRCh38 coordinates", file=sys.stderr)
            return 1
        variant = resolved
        print(f"Resolved {variant.rsid} -> {variant.hgvs_g}")

    # Step 1: VEP (required, resolves coordinates for HGVS coding inputs)
    vep, variant = query_vep(variant, proxy=args.proxy)
    if vep.get("error"):
        print(f"VEP query failed: {vep['error']}", file=sys.stderr)
        if not variant.chrom:
            return 1

    # Extract gene symbol for OMIM/UniProt
    gene_symbol = vep.get("gene_symbol")

    # Step 2: optional queries
    gnomad = {} if args.no_gnomad else query_gnomad(variant, proxy=args.proxy)
    clinvar = {} if args.no_clinvar else query_clinvar(variant, vep, proxy=args.proxy)
    omim = None if args.no_omim else query_omim(vep, omim_db=args.omim_db)
    uniprot = {}
    if gene_symbol:
        uniprot = query_uniprot(gene_symbol, protein_pos=vep.get("protein_start"), proxy=args.proxy)
    literature = {} if args.no_literature else query_literature(variant, vep, proxy=args.proxy)

    interpretation = build_interpretation(vep, clinvar)

    result = {
        "variant": {
            "raw": variant.raw,
            "chrom": variant.chrom,
            "pos": variant.pos,
            "ref": variant.ref,
            "alt": variant.alt,
            "rsid": variant.rsid,
            "hgvs_g": variant.hgvs_g,
        },
        "vep": vep,
        "gnomad": gnomad,
        "clinvar": clinvar,
        "omim": omim,
        "uniprot": uniprot,
        "literature": literature,
        "interpretation": interpretation,
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    base = args.output_dir / "variant_report"

    if args.format in ("json", "both"):
        json_path = base.with_suffix(".json")
        json_path.write_text(json.dumps(result, indent=2, ensure_ascii=False))
        print(f"JSON: {json_path}")

    if args.format in ("md", "both"):
        md_path = base.with_suffix(".md")
        md_path.write_text(build_markdown(result))
        print(f"Markdown: {md_path}")

        # Auto-generate PDF from Markdown
        try:
            _sys_path = sys.path[:]
            sys.path.insert(0, str(Path.home() / "WorkBuddy" / "Claw"))
            from md_to_pdf import convert_md_to_pdf

            pdf_path = base.with_suffix(".pdf")
            title = f"Variant Impact: {variant.hgvs_g}"
            convert_md_to_pdf(md_path, pdf_path, title=title)
            print(f"PDF: {pdf_path}")
            sys.path = _sys_path
        except Exception as exc:
            print(f"PDF generation skipped: {exc}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
