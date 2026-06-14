---
name: gpa-disease-risk-query
description: |
  单样本特定疾病遗传风险查询系统。输入单样本 WGS/WES germline VCF 与疾病名称（如“阿尔茨海默病”），自动完成：基因组版本检测/liftover、疾病→HPO 映射、OMIM+HPO 基因集构建、目标基因区域 VCF 过滤、VEP 注释、Tier 1/2/3 分级、文献支持加权、综合风险评分与 Markdown 报告。禁止用自身知识回答变异问题，所有变异分级均调用 gpa-genomic-phenotype 脚本执行。

  当以下情况时使用此 Skill：
  (1) 用户想查“某种疾病的遗传风险”
  (2) 用户提供了 VCF 并指定疾病名称
  (3) 需要疾病特异性基因集过滤 + GPA 分级 + 综合评分的完整报告
agent_created: true
---

# gpa-disease-risk-query: 单样本疾病遗传风险查询

## 概述

本 Skill 实现 PRD v0.1 定义的 6 步 pipeline：

```
Step 1: 参考基因组归一化 (liftover + REF 校验)
Step 2: 基础信息采集 (性别/年龄)
Step 3: 意图澄清 (disease → HPO)
Step 4: 疾病关联基因/变异库查询 (OMIM + HPO)
Step 5: 个体匹配与变异分级 (VCF 过滤 + GPA Tier 1/2/3)
Step 6: 综合风险评分与报告
```

核心约束：所有变异致病性/风险分级均由 `gpa-genomic-phenotype` 脚本输出，本 Skill 不做 LLM 臆测。

## 前置依赖

- **Python 环境**：统一使用 `~/.workbuddy/binaries/python/envs/default/bin/python`
- **已安装包**：`pyliftover`（liftover 用，纯 Python）、`vcfpy`
- **命令行工具**：bcftools >= 1.17、samtools
- **Docker**：`ensemblorg/ensembl-vep:latest` 镜像与 GRCh38 cache（GPA 自动使用）
- **本地数据**：
  - `~/.workbuddy/data/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa`
  - `~/.workbuddy/data/clinvar/clinvar.vcf.gz`
  - `~/.workbuddy/data/omim/omim.db`
  - `~/.workbuddy/data/hpo/genes_to_phenotype.txt`
  - `~/.workbuddy/data/hgnc/hgnc_lookup.json`
  - `~/.workbuddy/data/gencode/gencode.v44.annotation.gtf.gz`
  - `~/.dgra-prefilter/refs/gencode_v44_gene_loci.bed`
  - `~/.dgra-prefilter/refs/clinvar_pathogenic_GRCh38.bed`

## 快速开始

### 前置检查

```bash
python3 ~/.workbuddy/skills/gpa-disease-risk-query/scripts/main.py --preflight
```

### 运行完整查询

```bash
python3 ~/.workbuddy/skills/gpa-disease-risk-query/scripts/main.py \
  --vcf /path/to/sample.vcf.gz \
  --disease "Alzheimer disease" \
  --sex female \
  --age 45 \
  --output-dir /tmp/drq_output
```

### 指定 HPO ID

```bash
python3 ~/.workbuddy/skills/gpa-disease-risk-query/scripts/main.py \
  --vcf /path/to/sample.vcf.gz \
  --disease "阿尔茨海默病" \
  --hpo-id HP:0000726 \
  --sex female \
  --age 45
```

### 离线模式

```bash
python3 ~/.workbuddy/skills/gpa-disease-risk-query/scripts/main.py \
  --vcf /path/to/sample.vcf.gz \
  --disease "cardiomyopathy" \
  --offline
```

## CLI 参数

| 参数 | 必需 | 默认值 | 说明 |
|------|------|--------|------|
| `--vcf` / `--input` | 是 | — | 输入 VCF/VCF.gz |
| `--disease` | 是 | — | 疾病名称或自然语言描述 |
| `--hpo-id` | 否 | — | 显式 HPO ID |
| `--sex` | 否 | unknown | male / female / unknown |
| `--age` | 否 | — | 年龄 |
| `--family-history` | 否 | false | 一级亲属患病 |
| `--tissue` | 否 | auto | 覆盖 GPA 组织类型 |
| `--output-dir` | 否 | `./` | 输出目录 |
| `--max-genes` | 否 | 200 | 最大疾病基因数 |
| `--offline` | 否 | false | GPA 离线模式 |
| `--no-spliceai` | 否 | false | 关闭 SpliceAI |
| `--two-phase` | 否 | false | GPA 两阶段管线 |
| `--skip-liftover` | 否 | false | 即使检测到 GRCh37 也不 liftover |
| `--chain-file` | 否 | auto | 自定义 liftOver chain |
| `--literature-genes` | 否 | — | 逗号分隔文献基因 |
| `--literature-variants` | 否 | — | 文献变异 JSON 文件 |
| `--disease-mode` | 否 | auto | `mendelian` / `complex` / `auto`；见下方评分模式 |
| `--preflight` | 否 | — | 只运行依赖检查 |

## 输出

运行后在 `--output-dir` 下创建以疾病名命名的子目录，包含：

- `report.md` — 人类可读综合风险报告
- `result.json` — 结构化结果
- `normalized.vcf.gz` — 归一化后的 GRCh38 VCF
- `filtered_disease_genes.vcf.gz` — 目标基因区域子集 VCF
- `gpa_progress.jsonl` — GPA 进度日志

## 综合风险评分规则

| 维度 | 权重 | 规则 |
|------|------|------|
| Tier 1 | 40% | 每命中 +40；纯合截断/已知致病额外 +20 |
| Tier 2 | 25% | 每命中 +15；同基因 multi-hit 额外 +10 |
| 文献支持 | 20% | 位点级 +15；基因级 +5 |
| 人群罕见度 | 10% | gnomAD AF < 0.1% +10；0.1-1% +5；>1% 0 |
| 性别/年龄 | 5% | X 连锁隐性且男性 +5 |

| 总分 | 风险等级 |
|------|----------|
| ≥ 80 | 高风险 |
| 50–79 | 中等风险 |
| 20–49 | 低风险 |
| < 20 | 无明确风险 |

## 复杂表型：遗传贡献度评分模式

当 `--disease-mode complex`（或 `auto` 模式下内置模板/关键词判定为复杂疾病）时，评分逻辑从"致病性风险"切换为"遗传贡献度"，更适合高尿酸、痛风、糖尿病、肥胖、血脂异常、高血压等复杂表型。

| 维度 | 权重 | 规则 |
|------|------|------|
| 高外显率单基因变异 | 30% | Tier 1 每个 +30；Tier 2 每个 +15 |
| 罕见功能变异 | 25% | 核心基因内 rare missense +5、LoF +10；ClinVar LP/P +5；纯合 +3；AF 越低额外加分；单变异上限 20 |
| GWAS 风险等位基因 | 25% | 直接命中 lead SNP：杂合 +3、纯合 +6；该维度上限 25 |
| 文献/通路支持 | 10% | 变异级 +8、基因级 +3 |
| 人群罕见度 | 10% | 最低 gnomAD AF < 0.01% 满分；< 0.1% 7 分；< 1% 3 分 |

| 总分 | 贡献度等级 | 含义 |
|------|------------|------|
| ≥ 70 | 遗传贡献较高 | 遗传因素在表型中占比较重 |
| 40–69 | 遗传贡献中等 | 遗传是多个因素之一 |
| 20–39 | 遗传贡献较低 | 少量相关变异，不足以单独解释 |
| < 20 | 无明确遗传贡献 | 未检出实质遗传证据 |

报告会额外输出"关键突变/位点贡献明细"，列出每个罕见功能变异和 GWAS lead SNP 的贡献分与依据。

## 关于 LiftOver 的设计选择

**实现方案：pyliftover + samtools faidx + vcfpy。**

理由：
1. **优先使用现有工具**：`pyliftover` 是纯 Python 库，直接读取 UCSC chain 文件；`vcfpy` 读取/写入 VCF；`samtools faidx` 查询目标参考序列。
2. **CrossMap 在当前环境不可用**：default venv 中 `CrossMap` 的编译依赖 `pyBigWig` 因 macOS 代码签名/Team-ID 冲突无法加载，故未采用。
3. **chain 文件不适合维护成数据表**：chain 是 UCSC 区间映射，本质是大量区间对；转成关系表会冗余且查询慢。改为缓存 UCSC chain 文件到 `~/.workbuddy/data/liftover/`，首次自动下载，后续复用。用户也可通过 `--chain-file` 自定义。
4. **REF 校验与 allele 修正**：liftover 后使用目标 FASTA 重新校验 REF；对于 SNP，若目标参考碱基与原 ALT 一致，则自动交换 REF/ALT 并反转基因型编码。
5. **备选**：若未来需要更高吞吐，可启用 `bcftools +liftover` 插件或修复 CrossMap 的签名问题后切回。

## GPA VEP 115 兼容性处理

本地 `ensemblorg/ensembl-vep:latest` 为 VEP 115，已移除旧版 `--af_gnomad_exome` / `--af_gnomad_genome` 参数，新版参数为 `--af_gnomade` / `--af_gnomadg`，输出字段为 `gnomADe_AF` / `gnomADg_AF`。`gpa-genomic-phenotype` v0.10.5 仍使用旧参数，会导致 VEP 报错。

本 Skill 的兼容层：
1. **预注释**：在调用 GPA 前，先用 VEP 115 兼容参数对过滤后的子集 VCF 进行 Docker VEP 注释。
2. **Parser 补丁**：运行时 monkey-patch `dgra_input_parsers.VCFParser._csq_to_variant`，当传统 `gnomAD_AF` 不存在时，自动回退读取 `gnomADe_AF` / `gnomADg_AF`。
3. **GPA 跳过自带注释器**：传入预注释 VCF 后，GPA 检测到 `INFO/CSQ`，直接走解析分支，不再触发不兼容的 VEP 命令。

## 与现有技能的集成

| 现有技能 | 集成方式 |
|----------|----------|
| `dgra-prefilter` | 复用其预构建的 `gencode_v44_gene_loci.bed` 与 `clinvar_pathogenic_GRCh38.bed`；本 Skill 在疾病基因集上构建自定义 BED 并调用 `bcftools view -T` 硬过滤 |
| `gpa-genomic-phenotype` | 通过 `dgra_cli_wrapper.run_gpa_from_file()` 调用，但先由 `gpa_compat.py` 完成 VEP 115 预注释与 parser patch，再传入注释后的 VCF 进行 Tier 分级 |

## 异常处理

- VCF header 缺失 assembly：通过 chr1 长度推断；仍不确定则提示用户。
- LiftOver 后 REF 不匹配率 > 5%：中断并提示检查 VCF。
- HPO 无匹配：降级为 OMIM 关键词搜索，仍无结果则输出空报告。
- GPA 返回空 Tier 1/2：综合评分落入“无明确风险”，仍输出完整报告。
- 本地数据缺失：`--preflight` 会标红，必要时切换 `--offline`。

## 分析定位：基因型-表型关联，而非分子机制解释

本 Skill 的设计目标是 **genotype-phenotype association / 遗传风险贡献度评估**，不是解释“某个突变如何通过具体蛋白通路导致表型”。

- **我们做的**：把样本基因型与已报道的疾病/表型关联证据（OMIM、HPO、ClinVar、GWAS、文献）进行匹配，按证据强度加权，给出“遗传上有多大可能倾向该表型”。
- **我们不做的**：不做因果推断、不重建信号通路、不解释变异对蛋白结构或细胞功能的直接影响。

因此，报告中的“贡献分”反映的是**已知证据下的关联强度**，而不是分子层面的致病确定性。最终临床解读权在用户手中。

## 输入 VCF 约定

- 本 Skill 默认用户提供的 VCF 是 **genotyping 后的样本级 VCF**（即包含 reference、heterozygous、homozygous 三种状态的完整位点集合）。
- **未在 VCF 中检出的位点，默认视为参考基因型（REF/REF）**，而不是"数据缺失"。
- 若 VCF 在上游被过滤（例如 GVCF genotype 后去除了常见变异、或仅保留 rare variants），则某些重要位点（如 APOE rs7412/rs429358）可能完全不出现在 VCF 中。此时 pipeline 会在报告中明确标注该位点"未检出/无法推断"，而不是将其直接判为参考基因型。
- 建议：进行疾病风险查询前，尽量使用未经过滤的原始 germline VCF；如需评估晚发性 AD 等常见变异驱动的疾病，APOE 等关键位点必须保留在 VCF 中。

## 数据局限性

- 仅分析 germline SNV/Indel，不覆盖 CNV/SV/表观遗传。
- 文献模块 MVP 仅支持通过 `--literature-genes`/`--literature-variants` 传入结果；PubMed 自动检索待 Phase 2 扩展。
- 综合风险评分仅供科研参考，不构成临床诊断。

## 文件结构

```
gpa-disease-risk-query/
  SKILL.md                 # 本文件
  config.json              # Skill 元数据
  scripts/
    main.py                # CLI 入口
    pipeline.py            # 6 步 pipeline 编排
    liftover.py            # pyliftover liftover + REF 校验
    hpo_mapper.py          # disease → HPO 映射
    gene_set_builder.py    # OMIM + HPO 基因集
    vcf_filter.py          # 目标基因 BED 构建 + VCF 过滤
    gpa_runner.py          # 调用 gpa-genomic-phenotype
    gpa_compat.py          # VEP 115 兼容层（预注释 + parser patch）
    risk_scorer.py         # 综合风险评分
    report.py              # Markdown 报告
    constants.py           # 路径与参数常量
  references/              # 运行时缓存（gene coords, HPO map）
  tests/                   # 测试
```
