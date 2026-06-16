# GRCh38 Variant Functional Impact Analyzer

> 一键分析 GRCh38 单变异的功能影响 — 从 VEP 预测到 ClinVar 分类，从 gnomAD 频率到 OMIM 表型，自动生成 JSON / Markdown / PDF 三格式报告。

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

---

## 功能概览

输入一个变异标识（坐标 / HGVS / rsID），自动查询 6 个公共知识库，生成结构化注释报告：

| 数据源 | 查询内容 | 离线？ |
|--------|----------|--------|
| **Ensembl VEP** | SIFT, PolyPhen, AlphaMissense, CADD, SpliceAI, gnomAD 频率 | 需联网 |
| **gnomAD v4.1** | 外显子组 / 全基因组 AC/AN/AF | 需联网 |
| **NCBI ClinVar** | 临床意义分类, 审核状态, 关联表型 | 需联网 |
| **OMIM** | 基因-表型关联 (SQLite 本地查询) | ✅ 本地 |
| **UniProt** | 蛋白功能, 结构域, 变异附近特征 | 需联网 |
| **Europe PMC** | 变异特异性文献 | 需联网 |

---

## 输出

默认 `--format both` 生成三个文件：

| 文件 | 说明 |
|------|------|
| `variant_report.json` | 结构化机器可读结果 |
| `variant_report.md` | 中文 Markdown 综合报告 |
| `variant_report.pdf` | 可打印 PDF（含中文字体） |

---

## 前置依赖

### Python 环境

需要 Python ≥ 3.10，安装以下包：

```bash
pip install requests fpdf2
```

可选 — 如果需要从 rsID 解析坐标（Ensembl API 默认已包含）：

```bash
pip install myvariant   # 可选: rsID 解析备用方案
```

### 外部数据（可选）

**OMIM 本地数据库**：将 `omim.db` 放到 `~/.workbuddy/data/omim/omim.db`，或通过 `--omim-db` 指定路径。没有时 OMIM 章节自动跳过。

### 网络要求

需要访问以下境外 API（大陆用户请确保代理可用）：
- `rest.ensembl.org` — VEP + rsID 解析
- `gnomad.broadinstitute.org` — gnomAD GraphQL
- `eutils.ncbi.nlm.nih.gov` — ClinVar E-utilities
- `rest.uniprot.org` — UniProt
- `www.ebi.ac.uk` — Europe PMC

可使用 `--proxy http://127.0.0.1:7890` 指定 HTTPS 代理。

---

## 安装到 AI 编程助手

### WorkBuddy

**方法 1 — 从 Skill 市场安装**

在 WorkBuddy 对话中输入：

```
安装 grch38-variant-impact skill
```

**方法 2 — 手动安装**

```bash
git clone https://github.com/lzr098/Disease-Risk-Query.git \
  ~/.workbuddy/skills/grch38-variant-impact
```

然后在对话中通过 `@skill:grch38-variant-impact` 或自然语言触发（如"分析这个变异 chr2:219425737:G:C"）。

---

### Claude Code (Anthropic)

**方法 1 — 通过 CLAUDE.md 注册自定义命令**

在项目根目录的 `CLAUDE.md` 中添加：

```markdown
## 自定义技能: 变异功能分析

当用户要求分析某个变异的功能影响时，运行：

```bash
~/.workbuddy/skills/grch38-variant-impact/scripts/analyze_variant.py \
  --variant "{variant}" \
  --omim-db ~/.workbuddy/data/omim/omim.db \
  --output-dir ./variant_reports
```
```

Claude 会在对话中自动调用该脚本。

**方法 2 — 直接命令行调用**

在 Claude Code 会话中使用 Bash 工具执行：

```bash
python ~/.workbuddy/skills/grch38-variant-impact/scripts/analyze_variant.py \
  --variant "chr2:219425737:G:C" \
  --output-dir ./reports
```

---

### OpenAI Codex CLI

**方法 1 — 注册为 Codex Skill**

在 `~/.codex/skills/` 下创建软链接：

```bash
mkdir -p ~/.codex/skills
ln -s ~/.workbuddy/skills/grch38-variant-impact ~/.codex/skills/grch38-variant-impact
```

Codex 会自动发现 `SKILL.md` 中的 YAML frontmatter 并注册该技能。

**方法 2 — Shell 别名**

在 `~/.zshrc` 或 `~/.bashrc` 中添加：

```bash
alias variant-impact='python ~/.workbuddy/skills/grch38-variant-impact/scripts/analyze_variant.py'
```

---

### 其他支持 MCP 或自定义命令的 AI 工具

将脚本路径注册为工具即可。接受 `--variant`、`--output-dir` 等标准参数。

---

## 使用示例

### 命令行

```bash
python scripts/analyze_variant.py \
  --variant "chr2:219425737:G:C" \
  --omim-db ~/.workbuddy/data/omim/omim.db \
  --output-dir ./reports/des_p455h
```

### 支持的输入格式

```bash
# 坐标格式
--variant "chr2:219425737:G:C"
--variant "2:219425737:G:C"

# HGVS 基因组
--variant "2:g.219425737G>C"

# dbSNP rsID（自动解析到 GRCh38）
--variant "rs267607488"

# HGVS 编码（VEP 自动解析）
--variant "NM_001927.4:c.1363G>C"
```

### 可选参数

```
--format {json,md,both}   输出格式（默认 both）
--omim-db PATH             OMIM SQLite 数据库路径
--output-dir PATH          输出目录
--proxy URL                HTTPS 代理
--no-gnomad               跳过 gnomAD 查询
--no-clinvar              跳过 ClinVar 查询
--no-literature           跳过文献检索
--no-omim                 跳过 OMIM 查询
```

---

## 输出字段 (JSON)

```json
{
  "variant": {
    "raw": "原始输入",
    "chrom": "2", "pos": 219425737,
    "ref": "G", "alt": "C",
    "rsid": "rs267607488",
    "hgvs_g": "2:g.219425737G>C"
  },
  "vep": {
    "transcript": "ENST00000373960",
    "cdna": "c.1363G>C",
    "protein": "p.Asp455His",
    "sift": {"prediction": "deleterious", "score": 0},
    "polyphen": {"prediction": "probably_damaging", "score": 0.967},
    "alphamissense": {"class": "likely_pathogenic", "pathogenicity": 0.9102},
    "cadd_phred": 32,
    "spliceai": {"DS_AG": 0, "DS_AL": 0, "DS_DG": 0, "DS_DL": 0}
  },
  "gnomad": {
    "exome": {"ac": 0, "an": 1461760, "af": 0.0},
    "genome": {"ac": 0, "an": 152172, "af": 0.0}
  },
  "clinvar": {
    "accession": "VCV0000...",
    "classification": "Uncertain significance",
    "review_status": "criteria provided, ..."
  },
  "omim": {
    "gene_symbol": "DES",
    "phenotypes": [
      {"phenotype": "Cardiomyopathy, dilated, 1I", "inheritance": "AD"}
    ]
  },
  "uniprot": {
    "accession": "P17661",
    "protein_name": "Desmin",
    "features_near_variant": [...]
  },
  "literature": {
    "query": "...",
    "count": 12,
    "articles": [...]
  },
  "interpretation": {
    "summary": "VUS with suspicious in-silico evidence",
    "evidence": ["AlphaMissense: likely_pathogenic", "CADD phred 32", ...]
  }
}
```

---

## 限制

- **仅支持 GRCh38**，不包含 liftover 功能。
- **仅支持 SNV 和小 indel**，不支持 SV/CNV。
- **API 依赖网络**，境外 API 在大陆可能需要代理。
- **OMIM 需要本地 SQLite**，未提供时 OMIM 章节自动跳过。
- **解读仅供参考**，不可直接作为临床诊断依据。

---

## 疾病风险查询 — 内置模板管理

除了单变异分析，本仓库也提供**疾病级基因组风险查询** (`gpa-disease-risk-query` skill)。它基于预置疾病模板 + 多源基因集自动构建，完成从 VCF 到分层风险报告的完整流程。

### 已有内置模板

| 疾病 | 模式 | 基因数 | 已知致病位点 | GWAS SNP | 文献 |
|------|------|--------|-------------|----------|------|
| Alzheimer disease | mixed | 245 | 13 | 77 | 10 |
| Parkinson disease | mixed | 170 | 10 | 39 | 6 |
| Adult vision disorders | complex | 83 | 10 | 56 | 5 |
| Hyperuricemia / Gout | complex | 69 | 8 | 40 | 10 |
| Myocardial infarction / CAD | complex | 54 | 8 | 37 | 4 |

每个模板包含：
- **基因集** — 按 Mendelian / GWAS 分层，标注外显率、贡献权重、证据类型
- **已知致病位点** — 精选 ClinVar P/LP + OMIM + 文献位点，附基因型贡献说明
- **GWAS lead SNP** — 来自最大规模 GWAS 的效应等位基因、OR/beta、人群频率
- **关键蛋白域** — 突变热点结构注释（如 APOB 受体结合域、DES Tail 域）
- **关键文献** — PMID + 疾病关联证据

### 如何添加新疾病模板

将模板定义写入 `scripts/constants.py` 的 `DISEASE_BUILTIN_REFS` 字典。格式参考已有模板：

```python
DISEASE_BUILTIN_REFS = {
    ...
    "colon cancer": {
        "aliases": ["colon cancer", "colorectal cancer", "crc"],
        "mode": "complex",
        "gene_set": [
            {"gene": "APC", "tier": "mendelian_high", "contribution_score": 1.0,
             "penetrance": ">0.95", "penetrance_score": 0.95, "evidence": "familial",
             "note": "FAP / Lynch syndrome core gene"},
            {"gene": "TP53", "tier": "mendelian_high", ...},
            ...
        ],
        "known_pathogenic_variants": [...],
        "gwas_lead_snps": [...],
        "key_literature": [...],
    }
}
```

同时更新两个映射表：

1. **`DISEASE_NAME_ALIASES`** — 中文/英文别名 → 规范名
2. **`COMPLEX_DISEASE_KEYWORDS`** — 用于 `complex` 模式的疾病关键词列表

### 未匹配疾病的处理流程

当用户查询的疾病**没有**内置模板时，系统自动走动态构建路径：

```
用户输入 "亨廷顿病"
    │
    ├─ 1. 别名匹配 ─── 找到 hpo → 映射 HP:0002072
    │
    ├─ 2. HPO 基因提取 ─── genes_for_hpo("HP:0002072") → HTT, JPH3, ...
    │
    ├─ 3. OMIM 关键词检索 ─── 本地 omim.db 全文检索 "huntington"
    │
    ├─ 4. ClinVar 安全网 ─── 提取 P/LP 位点补充基因集
    │
    ├─ 5. 文献动态检索 ─── PubMed API 提取近年关键文献
    │
    └─ 6. DiseaseProfile 构建 → 自动评分 → 报告
```

**关键决策点：**

| 场景 | 行为 |
|------|------|
| 疾病名匹配内置模板 | 使用预置基因集，**不**额外查询 OMIM/HPO（保持模板整洁和快速） |
| 疾病名未匹配但有 HPO 映射 | 动态构建：HPO → OMIM → ClinVar → 文献 |
| 疾病名完全没有匹配 | 以用户输入作为关键词，直接从 OMIM + PubMed 检索 |
| 想要强制动态扩展 | 使用 `--enrich` 参数，即便是内置模板也会叠加 HPO/OMIM 基因 |

**常见中文输入映射：**

```
阿尔茨海默病 / AD       → alzheimer disease
帕金森 / PD             → parkinson disease
高尿酸血症 / 痛风        → hyperuricemia
心肌梗死 / 冠心病 / CAD  → myocardial infarction
视力下降 / 黄斑变性      → adult vision disorders
```

### 审查与迭代

模板通过独立的 **`validate_disease_templates.py`** 脚本进行自动化校验：

- 必需字段完整性检查（gene_set, known_pathogenic_variants 等）
- 重复基因/rsID 检测（允许同 rsID 不同等位基因）
- REF 等位基因一致性（vs GRCh38 FASTA）
- 贡献评分合理性范围检查

每次模板修改后运行验证：

```bash
cd ~/.workbuddy/skills/gpa-disease-risk-query
python scripts/validate_disease_templates.py
```

---

## 目录结构

```
grch38-variant-impact/
├── SKILL.md                   # 技能定义（YAML frontmatter）
├── README.md                  # 本文件
├── scripts/
│   └── analyze_variant.py     # 主分析脚本
└── references/
    └── api_reference.md       # API 端点参考文档
```

---

## 许可

MIT License. 详见 [LICENSE](LICENSE) 文件。
