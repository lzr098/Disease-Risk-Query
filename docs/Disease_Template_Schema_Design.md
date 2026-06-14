# gpa-disease-risk-query 疾病模板新 Schema 设计

## 目标

让疾病内置模板同时携带：
1. **基因层级**：对疾病的贡献度 / 外显率 / 证据类型；
2. **SNP 层级**：效应等位基因、beta/OR、EAF、置信度、贡献度；
3. **区域层级**：重要调控/编码区域；
4. **评分参与**：这些字段直接参与最终综合风险评分。

---

## 通用 Schema

```python
DISEASE_BUILTIN_REFS: dict[str, dict] = {
    "<disease_key>": {
        "mode": "mendelian" | "complex" | "auto",
        "gene_set": [
            {
                "gene": str,               # HGNC symbol
                "tier": str,               # mendelian_high / mendelian_mod / strong_gwas / gwas / supporting
                "contribution_score": float, # 0.0–1.0，该基因对疾病风险的相对贡献
                "penetrance": str,         # 人类可读：">0.95", "~0.3", "dosage", "very_low"
                "penetrance_score": float, # 0.0–1.0，用于量化计算
                "evidence": str,           # familial / gwas / rare_variant / functional
                "note": str,
            }
        ],
        "gwas_lead_snps": [
            {
                "rsid": str,
                "chrom": str,              # chrN，GRCh38
                "pos": int,                # 1-based
                "gene": str,
                "effect_allele": str,      # 风险/效应等位基因
                "other_allele": str,       # 对照等位基因
                "beta": float | None,      # log(OR)
                "or": float | None,        # odds ratio
                "eaf_eur": float | None,   # 欧洲效应等位基因频率
                "eaf_eas": float | None,   # 东亚效应等位基因频率
                "tier": str,               # high_effect / gwas / tag
                "contribution_score": float, # 0.0–1.0
                "confidence": str,         # high / moderate / low
                "note": str,
            }
        ],
        "regions": [
            {
                "chrom": str,
                "start": int,              # 0-based
                "end": int,
                "gene": str,
                "type": str,               # coding_exon / gwas_locus / regulatory / enhancer
                "contribution_score": float,
                "note": str,
            }
        ],
        "key_literature": [
            {
                "pmid": str,
                "title": str,
                "genes": list[str],
                "evidence_type": str,      # familial_penetrance / gwas / functional_review
                "note": str,
            }
        ],
    }
}
```

---

## 字段含义与使用方式

| 字段 | pipeline 中的用途 |
|---|---|
| `gene_set[].tier` | Tier 1/2/3 优先级：mendelian_high 直接进 Tier 1；gwas 仅作 Tier 3 降噪和注释。 |
| `gene_set[].contribution_score` | 落在该基因的变异额外加权；HIGH 影响变异 × 高贡献基因得分更高。 |
| `gene_set[].penetrance_score` | 用于报告和患者沟通，暂时不直接乘入评分（避免单样本误判）。 |
| `gwas_lead_snps[].effect_allele` | 校验 VCF 中风险等位基因方向，自动换算 dosage。 |
| `gwas_lead_snps[].beta` | 直接累加到多基因贡献分：score += beta × dosage。 |
| `gwas_lead_snps[].contribution_score` | 当 beta 缺失时，用该分数作为替代权重。 |
| `regions[]` | 对非编码变异（regulatory/enhancer/gwas_locus）进行注释和 proximity 打分。 |
| `key_literature[].evidence_type` | 文献支持打分时区分孟德尔外显率证据与 GWAS 证据。 |

---

## 兼容性策略

- `disease_reference.py` 在读取模板时做**归一化**：
  - 老模板（仅 `core_genes`/`gwas_loci_genes`）自动转成新 `gene_set`，tier 和 contribution 用默认值。
  - 新模板直接解析。
- `pipeline.py` 中优先使用 `gene_set` / `gwas_lead_snps` 结构化数据；老字段作为 fallback。

---

## AD 模板 contribution_score 设计原则

| tier | contribution_score 范围 | 说明 |
|---|---|---|
| `mendelian_high` | 0.90–1.00 | APP/PSEN1/PSEN2 致病变异几乎决定疾病。 |
| `mendelian_mod` | 0.50–0.80 | APOE、TREM2、SORL1 等强风险基因，但非完全外显。 |
| `strong_gwas` | 0.30–0.50 | BIN1、PICALM、CLU、CR1 等经多次 GWAS 验证。 |
| `gwas` | 0.10–0.30 | 其他 PGS004863 位点。 |
| `supporting` | 0.05–0.10 | 间接或新兴证据。 |

SNP contribution_score 按 `abs(beta)` 或 `abs(OR-1)` 归一化到 0.05–0.40。
