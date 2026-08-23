# GPA Disease Risk Query &nbsp;·&nbsp; 单样本疾病遗传风险查询

[![version](https://img.shields.io/badge/version-2.2.0-blue)](https://github.com/lzr098/Disease-Risk-Query)
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
- Supports **pre-built disease templates** for 19 diseases
- Gene-level **inheritance mode** (AD/AR/XL/MT/SD) from ClinGen + HPO for zygosity interpretation
- Dynamic gene set construction for unmatched diseases
- China-population reference (CHN_AF from local VCF)

## 做了什么

- 输入**疾病名称** (中文/英文) → HPO 映射 → OMIM + HPO + ClinVar 基因集构建
- 从 VCF 中提取目标基因区域变异
- VEP 功能注释 (SIFT, PolyPhen, CADD, SpliceAI, REVEL, gnomAD 频率)
- **Tier 1** (高风险) / **Tier 2** (中风险) / **Tier 3** (低风险) 三级分层
- ClinVar 致病/可能致病变异 + OMIM 表型 + GWAS lead SNP + 关键文献佐证
- 支持 **19 个内置疾病模板** (阿尔茨海默病、帕金森、心梗、2型糖尿病、乳腺癌、结直肠癌、IBD、骨质疏松、FH、MDS、AGA、前列腺癌、精神疾病、视力障碍等)
- 基因级**遗传模式** (AD/AR/XL/MT/SD，来自 ClinGen + HPO) 支持合子状态解读
- 未匹配疾病**动态基因集构建**
- 中国人群频率参考 (CHN_AF)

---

## Quick Start · 快速开始

```bash
python scripts/main.py \
  --vcf patient.vcf.gz \
  --disease "阿尔茨海默病" \
  --output-dir ./alzheimer_report
```

```bash
python scripts/main.py \
  --vcf patient.vcf.gz \
  --disease "breast cancer" \
  --output-dir ./breast_cancer_report \
  --refresh-cache
```

---

## Built-in Disease Templates · 内置疾病模板

共 **19 个**内置模板，`mendelian` 基因均带 `inheritance` 遗传模式字段（ClinGen + HPO 数据库）。

| Disease · 疾病 | Mode · 模式 | 基因数 | Mendelian¹ | GWAS |
|---|---|---|---|---|
| Alzheimer disease · 阿尔茨海默病 | mixed | 245 | 8 | 237 |
| Psychiatric disorders · 精神疾病 (BPD/MDD/ADHD/PTSD/Bipolar) | complex | 184 | 12 | 172 |
| Parkinson disease · 帕金森病 | mixed | 170 | 16 | 154 |
| Adult vision disorders · 成人视力障碍 | complex | 83 | 37 | 46 |
| Circadian glucose regulation · 昼夜血糖调节 | complex | 66 | 0 | 66 |
| Hyperuricemia / Gout · 高尿酸血症 | complex | 61 | 42 | 19 |
| FSHD complications · FSHD 并发症 | complex | 58 | 29 | 29 |
| Myocardial infarction / CAD · 心肌梗死 | complex | 52 | 11 | 41 |
| Body mass index · BMI | complex | 39 | 11 | 28 |
| Androgenetic alopecia · 雄激素性脱发 | complex | 37 | 1 | 36 |
| Colorectal cancer · 结直肠癌 | complex | 33 | 18 | 15 |
| Breast cancer · 乳腺癌 | complex | 32 | 16 | 16 |
| Familial hypercholesterolemia · 家族性高胆固醇血症 | complex | 32 | 9 | 23 |
| Osteoporosis · 骨质疏松 | complex | 32 | 10 | 22 |
| Type 2 diabetes · 2型糖尿病 | complex | 31 | 9 | 22 |
| IBD · 炎症性肠病 | complex | 29 | 7 | 22 |
| Myelodysplastic syndromes · 骨髓增生异常综合征 | mixed | 24 | 20 | 4 |
| Prostate cancer · 前列腺癌 | complex | 23 | 20 | 3 |
| Retinal vasculopathy · 视网膜血管病 | mendelian | 10 | 8 | 2 |

> ¹ Mendelian 层（`mendelian_high` + `mendelian_mod`）基因，均带 `inheritance` 遗传模式字段。
> 多遗传模式基因（如 BRCA1 的 AD+AR）用逗号连接，不硬猜主模式。

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
| **PRS High** | High-weight PRS variants (weight 0.9) for top GWAS loci |
| **Literature Evidence** | PubMed PMID count, publication recency |
| **Inheritance Mode** | AD / AR / XL / XLd / XLr / MT / YL / SD — from ClinGen gene-disease validity + HPO (`references/gene_inheritance.json`, 5176 基因) |

---

## Key Options · 关键参数

| Flag | Description |
|---|---|
| `--vcf PATH` | Input germline VCF (required) |
| `--disease "name"` | Disease name, Chinese or English (required) |
| `--output-dir PATH` | Output directory |
| `--offline` | Skip API queries, use only local data |
| `--refresh-cache` | Force rebuild of disease reference cache |
| `--disease-mode MODE` | Scoring mode: mendelian / complex / auto |
| `--sex male\|female` | Sample sex (affects X-linked scoring) |
| `--assume-genotyped` | Treat missing sites as ref/ref (genotyped VCF) |

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
| **Gene inheritance** | **Local** | **`references/gene_inheritance.json` (ClinGen + HPO 合并, 5176 基因)** |
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

## Changelog · 更新日志

### v2.2.0 (2026-08-23)

- **基因遗传模式 (inheritance) 字段**：19 个模板的 Mendelian 层基因全部补上 `inheritance` 字段（AD/AR/XL/XLd/XLr/MT/YL/SD），来自 ClinGen gene-disease validity + HPO 合并数据库
- **新增 `references/gene_inheritance.json`**：ClinGen + HPO 合并的遗传模式数据库（5176 基因，285KB）
- **新增 `scripts/backfill_inheritance.py`**：可复用的补字段脚本（支持 `--dry-run`），规则：仅 Mendelian 层、ClinGen 优先 HPO 兜底、排除 Disputed/Refuted、多值逗号连接不硬猜
- **评分器修复**：`_af_rarity_factor()` AF 三态处理（AF=None 降权 ×0.5 而非当"罕见"）+ `_clinvar_evidence_level()`（仅 P/LP + VUS 的 missense 计中等外显，VUS 降权）——修复了"成年视力下降 74 分虚高"问题
- **build 检测修复**：支持 b37/b38 标记，修正 chr1 长度 fallback 判定逻辑（原逻辑写反）
- **模板数**：14 → 19（新增 psychiatric disorders、circadian glucose regulation、body mass index、FSHD complications、retinal vasculopathy）

### v2.1.0 (2026-07-10)

- **New template**: Androgenetic alopecia (AGA) — 37 genes, 9 GWAS SNPs, 4 PRS high, 10 genomic regions (Pirastu 2017)
- **New template**: Prostate cancer — 23 genes, 6 GWAS SNPs, 7 ClinVar P/LP
- **Bug fix**: `variant_class` field missing in `prs_variants_high` across 5 templates (MI, Alzheimer, T2D, Hyperuricemia, AGA) — `VariantWeight.from_dict()` silently ignored `'tier'` key, causing prs_high entries (weight 0.9) to be scored as gwas_prs (weight 0.3)
- **AGA template optimization**: gene_set 16→37, added prs_variants_high (4), prs_variants (3), regions (10), core_genes (8), gwas_loci_genes (34), key_regions 3→8, key_literature 3→7, beta=ln(OR) for all GWAS SNPs
- Template count: 12 → 14

### v2.0.0

- Initial public release with 12 disease templates
- VEP Docker offline annotation support
- Two-phase GPA pipeline
- ReportComposer integration

---

## License

MIT

---

**Maintainer**: [@lzr098](https://github.com/lzr098)
