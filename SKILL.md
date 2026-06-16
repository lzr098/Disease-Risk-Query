---
name: grch38-variant-impact
description: Analyze the functional impact of a single GRCh38 variant by querying Ensembl VEP, gnomAD, ClinVar, local OMIM, UniProt, and Europe PMC. Returns a structured JSON/Markdown report with pathogenicity predictions, population frequencies, gene-phenotype associations, protein context, and literature evidence.
agent_created: true
---

# GRCh38 Variant Functional Impact Analyzer

## Overview

This skill performs a one-stop functional annotation for a single GRCh38 variant. It aggregates publicly available evidence from multiple sources and produces a structured report suitable for variant interpretation.

## When to Use

Use this skill when the user asks to:

- Analyze the functional impact of a specific variant.
- Look up ClinVar, gnomAD, VEP, or OMIM annotations for a GRCh38 variant.
- Get a quick variant interpretation summary (VUS / benign / pathogenic evidence).
- Understand the gene/protein context of a missense or regulatory variant.

## Supported Input Formats

Pass the variant via `--variant` as one of:

1. **Chrom:Pos:Ref:Alt** — `chr2:21007456:G:C` or `2:21007456:G:C`
2. **HGVS genomic** — `2:g.21007456G>C`
3. **dbSNP rsID** — `rs755753065`
4. **HGVS coding** — `NM_000384.3:c.9412C>G` (VEP will resolve to genomic)

## Workflow

1. Parse the input variant into a canonical representation.
2. Query **Ensembl VEP** for functional predictions (SIFT, PolyPhen, AlphaMissense, CADD, SpliceAI) and embedded gnomAD frequencies.
3. Query **gnomAD v4.1 GraphQL** for overall exome/genome allele counts and frequencies.
4. Query **NCBI ClinVar** (E-utilities) for classification, review status, associated phenotypes, and last evaluation date.
5. Query the **local OMIM SQLite database** (`~/.workbuddy/data/omim/omim.db`) for gene-phenotype associations.
6. Query **UniProt** for protein length, function description, and features near the affected residue.
7. Search **Europe PMC** for publications specifically mentioning the variant.
8. Output a JSON file and a human-readable Markdown summary.

## Main Script

Run the analyzer directly:

```bash
python ~/.workbuddy/skills/grch38-variant-impact/scripts/analyze_variant.py \
  --variant "chr2:21007456:G:C" \
  --omim-db ~/.workbuddy/data/omim/omim.db \
  --output-dir ./variant_report
```

Optional flags:

- `--format {json,md,both}` — default `both`. When Markdown is generated, a **PDF** is automatically produced.
- `--no-gnomad` — skip gnomAD direct query
- `--no-clinvar` — skip ClinVar query
- `--no-literature` — skip Europe PMC search

## Output Files

When `--format both` (default):

- `variant_report.json` — structured machine-readable result
- `variant_report.md` — human-readable Markdown summary
- `variant_report.pdf` — printable PDF report (auto-generated from Markdown)

## Output Schema (JSON)

```json
{
  "input": {"raw": "chr2:21007456:G:C", "chrom": "2", "pos": 21007456, "ref": "G", "alt": "C", "hgvs_g": "2:g.21007456G>C"},
  "vep": {
    "transcript": "ENST00000233242",
    "cdna": "c.9412C>G",
    "protein": "p.Pro3138Ala",
    "sift": {"prediction": "deleterious", "score": 0},
    "polyphen": {"prediction": "benign", "score": 0.143},
    "alphamissense": {"class": "likely_benign", "pathogenicity": 0.1959},
    "cadd_phred": 19.4,
    "spliceai": {"DS_AG": 0, "DS_AL": 0, "DS_DG": 0, "DS_DL": 0},
    "gnomad_frequencies": {"gnomade": 8.9e-06, "gnomadg": 1.3e-05, "gnomade_eas": 2.27e-04, "gnomadg_eas": 3.85e-04}
  },
  "gnomad": {"exome": {"ac": 13, "an": 1461760, "af": 8.9e-06}, "genome": {"ac": 2, "an": 152172, "af": 1.3e-05}},
  "clinvar": {
    "accession": "VCV001766890",
    "classification": "Likely benign",
    "review_status": "criteria provided, multiple submitters, no conflicts",
    "last_evaluated": "2024/09/06",
    "traits": ["Cardiovascular phenotype", "Familial hypobetalipoproteinemia 1", "Hypercholesterolemia, autosomal dominant, type B"]
  },
  "omim": {
    "mim_number": "107730",
    "title": "Apolipoprotein B-100; APOB",
    "phenotypes": [
      {"phenotype": "Hypercholesterolemia, familial, 2", "mim_number": "144010", "inheritance": "AD"},
      {"phenotype": "Hypobetalipoproteinemia", "mim_number": "615558", "inheritance": "AR"}
    ]
  },
  "uniprot": {"accession": "P04114", "length": 4563, "function": "...", "features_near": []},
  "literature": {"query": "...", "count": 0, "articles": []},
  "interpretation": {
    "summary": "Likely benign / VUS-low",
    "evidence": ["ClinVar: likely benign", "AlphaMissense: likely_benign", "CADD phred 19.4"]
  }
}
```

## Important Limitations

- **GRCh38 only**. The skill does not lift over from GRCh37.
- **SNVs and small indels only**. Large SVs/CNVs are not supported.
- **External API dependency**. VEP, gnomAD, ClinVar, UniProt, and Europe PMC require internet access.
- **OMIM requires the local SQLite database** at `~/.workbuddy/data/omim/omim.db`. If missing, the OMIM section is skipped.
- **Interpretation is not a clinical diagnosis**. The generated summary is for research/screening purposes only.

## Scripts

- `scripts/analyze_variant.py` — main analyzer

## References

- `references/api_reference.md` — summary of queried APIs and their response fields
