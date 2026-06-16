# API Reference for grch38-variant-impact

This document summarizes the external APIs and local database used by the skill.

## Ensembl REST API

### Resolve rsID to coordinates

```http
GET https://rest.ensembl.org/variation/human/{rsid}?pops=0
Content-Type: application/json
```

Response field of interest: `mappings[].{seq_region_name, start, allele_string}`.

### Variant Effect Predictor (VEP)

```http
POST https://rest.ensembl.org/vep/human/hgvs
Content-Type: application/json

{
  "hgvs_notations": ["2:g.21007456G>C"],
  "AlphaMissense": true,
  "CADD": true,
  "SpliceAI": true
}
```

Key response fields:

| Field | Meaning |
|---|---|
| `transcript_consequences[].gene_symbol` | HGNC gene symbol |
| `transcript_consequences[].transcript_id` | Ensembl transcript ID |
| `transcript_consequences[].hgvsc` / `hgvsp` | cDNA / protein HGVS |
| `transcript_consequences[].sift_prediction` / `sift_score` | SIFT prediction |
| `transcript_consequences[].polyphen_prediction` / `polyphen_score` | PolyPhen prediction |
| `transcript_consequences[].alphamissense.am_class` | AlphaMissense class |
| `transcript_consequences[].alphamissense.am_pathogenicity` | AlphaMissense score |
| `transcript_consequences[].cadd_phred` | CADD phred-like score |
| `transcript_consequences[].spliceai.DS_*` | SpliceAI delta scores |
| `colocated_variants[].frequencies.{alt}.{gnomade,gnomadg,...}` | gnomAD frequencies |

## gnomAD GraphQL API

```http
POST https://gnomad.broadinstitute.org/api/
Content-Type: application/json
```

Query type for dataset variable: `DatasetId!` (e.g., `gnomad_r4`).

Key response fields:

| Field | Meaning |
|---|---|
| `region.variants[].variant_id` | gnomAD variant ID (`chrom-pos-ref-alt`) |
| `region.variants[].exome.{ac,an,af}` | Exome allele count / number / frequency |
| `region.variants[].genome.{ac,an,af}` | Genome allele count / number / frequency |

## NCBI E-utilities (ClinVar)

Search:

```http
GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term={term}&retmode=json&retmax=1
```

Summary:

```http
GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={uid}&retmode=json
```

Key response fields:

| Field | Path |
|---|---|
| Classification | `germline_classification.description` |
| Review status | `germline_classification.review_status` |
| Last evaluated | `germline_classification.last_evaluated` |
| Traits | `germline_classification.trait_set[].trait_name` |

## OMIM Local SQLite

Default path: `~/.workbuddy/data/omim/omim.db`

Relevant columns:

| Column | Use |
|---|---|
| `mim_number` | OMIM entry ID |
| `title` | Entry title |
| `hgnc_approved_gene_symbol` | Approved gene symbol (note: may contain leading/trailing spaces; use `TRIM()`) |
| `symbols` | Alternate symbols |
| `geneMap` | JSON list of gene-phenotype associations |

## UniProt REST API

Search gene to accession:

```http
GET https://rest.uniprot.org/uniprotkb/search?query=gene:{symbol}+AND+organism_id:9606+AND+reviewed:true&format=json&size=1
```

Fetch entry:

```http
GET https://rest.uniprot.org/uniprotkb/{accession}.json
```

Key response fields:

| Field | Path |
|---|---|
| Function | `comments[].commentType == FUNCTION` |
| Features | `features[].{type, description, location.start.value, location.end.value}` |
| Protein length | `sequence.length` |

## Europe PMC

```http
GET https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={terms}&format=json&pageSize=10
```

Response fields:

| Field | Path |
|---|---|
| Total hits | `hitCount` |
| Articles | `resultList.result[].{title, authorString, journalTitle, pubYear, pmid, doi}` |
