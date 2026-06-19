# GPA Disease Risk Query &nbsp;·&nbsp; 单样本疾病遗传风险查询

[![version](https://img.shields.io/badge/version-2.0.0-blue)](https://github.com/lzr098/Disease-Risk-Query)
[![python](https://img.shields.io/badge/python-3.10%2B-green)](https://www.python.org/)
[![license](https://img.shields.io/badge/license-MIT-lightgrey)](LICENSE)

> **Single-sample disease-specific genetic risk analysis.** From a germline WGS/WES VCF + disease name, automatically builds gene sets, annotates variants, classifies by Tier 1/2/3, supplements with literature evidence, and generates structured reports.

> **单样本特定疾病遗传风险查询。** 输入 germline VCF + 疾病名称，自动构建基因集、注释变异、Tier 1/2/3 分级、文献佐证，输出结构化风险报告。

---

## What It Does

- Maps disease names to **HPO terms** → builds gene sets from OMIM + HPO + ClinVar
- Filters VCF to disease-relevant genomic regions
- Annotates with Ensembl VEP (SIFT, PolyPhen, CADD, SpliceAI, REVEL, gnomAD AF)
- Classifies variants into **Tier 1** (high risk) / **Tier 2** (moderate) / **Tier 3** (low)
- Enriches with ClinVar P/LP, OMIM phenotypes, GWAS lead SNPs, key literature
- Supports **pre-built disease templates** for 10+ common diseases
- Dynamic gene set construction for unmatched diseases
- China-population reference (CHN_AF from local VCF)

## 做了什么

- 输入**疾病名称** (中文/英文) → HPO 映射 → OMIM + HPO + ClinVar 基因集构建
- 从 VCF 中提取目标基因区域变异
- VEP 功能注释 (SIFT, PolyPhen, CADD, SpliceAI, REVEL, gnomAD 频率)
- **Tier 1** (高风险) / **Tier 2** (中风险) / **Tier 3** (低风险) 三级分层
- ClinVar 致病/可能致病变异 + OMIM 表型 + GWAS lead SNP + 关键文献佐证
- 支持 **10+ 内置疾病模板** (阿尔茨海默病、帕金森、心梗、2型糖尿病、乳腺癌、结直肠癌等)
- 未匹配疾病**动态基因集构建**
- 中国人群频率参考 (CHN_AF)

---

## Quick Start · 快速开始

```bash
python scripts/run_P010.py \
  --vcf patient.vcf.gz \
  --disease "阿尔茨海默病" \
  --output-dir ./alzheimer_report
```

```bash
python scripts/run_P010.py \
  --vcf patient.vcf.gz \
  --disease "breast cancer" \
  --output-dir ./breast_cancer_report \
  --proxy http://127.0.0.1:7890
```

---

## Built-in Disease Templates · 内置疾病模板

| Disease · 疾病 | Mode · 模式 | Genes · 基因数 | GWAS SNPs | ClinVar P/LP |
|---|---|---|---|---|
| Alzheimer disease · 阿尔茨海默病 | mixed | 245 | 77 | 13 |
| Parkinson disease · 帕金森病 | mixed | 170 | 39 | 10 |
| Adult vision disorders · 成人视力障碍 | complex | 83 | 56 | 10 |
| Myocardial infarction / CAD · 心肌梗死 | complex | 54 | 37 | 8 |
| Type 2 diabetes · 2型糖尿病 | complex | 49 | 63 | 6 |
| Hyperuricemia / Gout · 高尿酸血症 | complex | 69 | 40 | 8 |
| Breast cancer · 乳腺癌 | complex | 32 | 44 | 12 |
| Colorectal cancer · 结直肠癌 | complex | 33 | 36 | 18 |
| IBD · 炎症性肠病 | complex | 29 | 46 | 3 |
| Osteoporosis · 骨质疏松 | complex | 32 | 49 | 5 |
| Familial hypercholesterolemia · 家族性高胆固醇血症 | mendelian | 32 | 18 | 38 |

---

## Pipeline Stages · 管线流程

```
Input: VCF + Disease name
  │
  ├─ 1. Disease → HPO mapping
  │     OMIM + HPO + ClinVar → gene set
  │
  ├─ 2. Pre-built template lookup
  │     matched → use template gene set + GWAS SNPs + P/LP
  │     unmatched → dynamic OMIM/HPO/ClinVar construction
  │
  ├─ 3. VCF filtering (disease-relevant regions)
  │
  ├─ 4. VEP annotation
  │     SIFT / PolyPhen / CADD / SpliceAI / REVEL / gnomAD AF
  │
  ├─ 5. Scoring & Tier classification
  │     ClinVar pathogenicity × gnomAD frequency × functional impact × domain × expression
  │
  ├─ 6. Literature enrichment
  │     Europe PMC + known pathogenic variant check
  │
  └─ 7. Report generation
        JSON + Markdown + disease risk summary
```

---

## Gene Contribution Scoring · 基因贡献评分

Each gene is scored on multiple axes:

| Axis · 维度 | Factor · 因子 |
|---|---|
| **Pathogenic Association** | ClinVar P/LP density, OMIM phenotype match |
| **Key Domains** | Missense/frameshift in critical protein domains → bonus |
| **ClinGen Validity** | Definitive × 1.0, Strong × 0.85, Moderate × 0.70, Limited × 0.50 |
| **GWAS Signal** | Lead SNP effect size (OR/beta), population AF |
| **Literature Evidence** | PubMed PMID count, publication recency |
| **Inheritance Match** | Dominant vs recessive vs X-linked pattern |

---

## Key Options · 关键参数

| Flag | Description |
|---|---|
| `--vcf PATH` | Input germline VCF (required) |
| `--disease "name"` | Disease name, Chinese or English (required) |
| `--output-dir PATH` | Output directory |
| `--proxy URL` | HTTPS proxy for API access |
| `--offline` | Skip API queries, use only local data |
| `--enrich` | Force dynamic gene set expansion even for built-in templates |

---

## Data Sources · 数据来源

| Source | Type | Location |
|---|---|---|
| Ensembl VEP | API | REST + VEP Docker |
| gnomAD v4.1 | API | GraphQL |
| OMIM | Local | `~/.workbuddy/data/omim/omim.db` |
| ClinVar | Local | `~/.workbuddy/data/clinvar/clinvar.vcf.gz` |
| ClinGen | Local | `~/.workbuddy/data/clingen/` |
| HPO | API | Ontology lookup |
| Europe PMC | API | Literature search |
| GWAS Catalog | API | Lead SNP retrieval |

---

## Requirements · 运行环境

- Python 3.10+
- `requests`, `aiohttp`
- `~/.workbuddy/data/omim/omim.db` (OMIM SQLite, optional but recommended)
- `~/.workbuddy/data/clinvar/clinvar.vcf.gz` (local ClinVar VCF, optional)

---

## Related Skills · 相关技能

| Skill | Repo | Purpose |
|---|---|---|
| **grch38-variant-impact** | [lzr098/variant-impact](https://github.com/lzr098/variant-impact) | Single variant ACMG classification |
| **GPA** | [lzr098/dgra-genomic-risk](https://github.com/lzr098/dgra-genomic-risk) | Whole-genome phenotype association |
| **GPA Filter** | [lzr098/GPA-Filter](https://github.com/lzr098/GPA-Filter) | Genomic region pre-filter |
| **sensory-genomics** | [lzr098/sensory-genomics-skill](https://github.com/lzr098/sensory-genomics-skill) | Five-sense genetic analysis |

---

## License

MIT

---

**Maintainer**: [@lzr098](https://github.com/lzr098)
