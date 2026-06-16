# 疾病模板创建规范

> 本规范基于 2026-06-16 对 5 个内置模板的全面审查与优化制定。
> 所有新增疾病模板必须通过 `scripts/validate_disease_templates.py` 验证。

---

## 1. 最小模板结构

```python
DISEASE_BUILTIN_REFS = {
    "disease_name": {
        "aliases": ["canonical_name", "中文名", "缩写", ...],  # ≥3 个别名
        "mode": "complex",        # mendelian | complex | mixed

        "gene_set": [...],                     # 基因列表
        "known_pathogenic_variants": [...],    # KPV
        "gwas_lead_snps": [...],               # GWAS lead SNPs
        "prs_variants": [...],                 # 标准 PRS SNPs
        "prs_variants_high": [...],            # 高置信 PRS SNPs (≥3)

        "regulatory_regions": [],
        "key_regions": {...},                  # ≥2 基因的蛋白域注释
        "key_literature": [...],               # ≥3 篇 (≥1 篇 2020+)
    }
}
```

## 2. 基因分层规则

| tier | 条件 | CS 范围 | 示例 |
|------|------|--------|------|
| `mendelian_high` | 直接因果，外显率 >0.9 | 1.0 | LDLR(FH), APP(AD) |
| `mendelian_mod` | 直接因果，外显率 0.3-0.9 或剂量效应 | 0.6-0.85 | SORL1, LPA |
| `strong_gwas` | GWAS >100K 样本，跨族裔验证 | 0.22-0.30 | CDKN2B(9p21) |
| `gwas` | 标准 GWAS 位点 | 0.12-0.22 | 次要位点 |

**基因准入规则：**
- 只有**直接因果**（primary metabolic/signaling pathway）基因能进 mendelian tier
- 间接/继发性基因（如代谢病导致的高尿酸、脂营养不良导致的 CAD）→ 进 gwas 或删除
- 与疾病完全无关的基因（如干扰素病基因）→ 直接删除
- 删除前运行独立来源审计（HPO + OMIM + 文献）

## 3. 基因条目格式

每个基因必须包含以下字段：

```python
{
    "gene": "SYMBOL",
    "tier": "mendelian_high",
    "contribution_score": 1.0,
    "penetrance": ">0.95",
    "penetrance_score": 0.95,
    "evidence": "familial",
    "note": "功能简述",
    "phenotype_assoc": "疾病-基因关联 + 遗传模式",   # 必填
    "key_domains": "关键结构域(残基范围)",           # 必填
    "clingen_validity": "Definitive",               # 建议填
}
```

## 4. GWAS SNP 格式

每条 GWAS SNP 必须包含：

```python
{
    "rsid": "rs1333049",
    "chrom": "chr9",
    "pos": 22125504,
    "gene": "CDKN2B",
    "effect_allele": "C",       # 效应等位基因
    "other_allele": "G",
    "or_value": 1.29,           # 与 beta 至少填一个
    "beta": None,
    "eaf_eur": 0.47,            # 必填（欧洲频率）
    "eaf_eas": 0.52,            # 必填（东亚频率）
    "tier": "strong_gwas",
    "contribution_score": 0.22, # 按 tier 规则
    "confidence": "high",
    "note": "基因-位点-研究来源",
}
```

**获取 EAF**: 通过 Ensembl API `?pops=1` 查询，提取 `gnomADg:nfe` + `gnomADg:eas` 频率。

**跨疾病 CS 一致性**: 同一 tier 基因的 CS 标准差 ≤0.15。当前基准值见第 2 节。

## 5. 已知致病位点格式

```python
{
    "rsid": "rs5742904",
    "chrom": "chr2",
    "pos": 21006288,           # GRCh38 Ensembl 坐标 (非 dbSNP)
    "ref": "C",                # GRCh38 参考等位基因
    "alt": "T",
    "gene": "APOB",
    "effect_allele": "T",
    "other_allele": "C",
    "or_value": None,
    "beta": None,
    "variant_class": "known_pathogenic",
    "contribution_score": 1.0,
    "confidence": "high",
    "note": "蛋白变化; 遗传模式; 杂合/纯合贡献说明",
}
```

**校验流程**:
1. 用 rsID 查询 Ensembl API → 获取 GRCh38 坐标和 REF/ALT
2. 用 `samtools faidx` 验证 REF = GRCh38 FASTA 参考等位基因
3. 查询 ClinVar local VCF 或在线 E-utilities 验证 P/LP 分类
4. note 必须描述不同基因型的贡献差异

## 6. PRS 分层设计

| 层 | 字段 | 权重 | SNP 数 | 来源要求 |
|----|------|------|--------|---------|
| 高置信 | `prs_variants_high` | 0.9 | 3-5 | 已验证 PRS (PMID+独立队列) |
| 标准 | `prs_variants` | 0.3 | 20-50 | GWAS 汇总统计 top 位点 |

## 7. key_regions 格式

至少 2 个基因，每个基因需标注：

```python
"GENE": {
    "note": "基因功能描述",
    "regions": [
        {"name": "Domain", "residues": "start-end", "note": "域功能"},
    ],
    "critical_residues": [
        {"residue": "AminoAcid+Pos", "note": "突变影响"},
    ],
}
```

## 8. 关键文献格式

≥3 篇，至少 1 篇 2020 年后：

```python
{
    "pmid": "12345678",
    "title": "...",
    "genes": ["GENE1", "GENE2"],
    "note": "期刊+年份+要点",
    "evidence_type": "gwas|experimental|review",
}
```

## 9. 独立来源审计（新模板前必做）

1. 运行 HPO 基因提取 → 发现模板未覆盖的基因
2. 运行 OMIM 关键词检索 → 发现文献提及的基因
3. 对比 ClinVar P/LP 安全网 → 发现临床已知致病变异
4. 审核 OMIM 独有基因：是否为噪音匹配？

## 10. 验证步骤

```bash
# 1. 语法检查
python -m py_compile scripts/constants.py

# 2. 模板验证
python scripts/validate_disease_templates.py

# 3. 全量回归
pytest tests/ -q
```
