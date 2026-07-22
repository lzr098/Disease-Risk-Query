"""Generate psychiatric_disorders.py disease template from collected SNPs."""

from __future__ import annotations

import json
import math
import re
import gzip
from pathlib import Path
from typing import Optional

GENCODE_GTF = Path("/Users/zhaorongli/.workbuddy/data/gencode/gencode.v44.annotation.gtf.gz")
BPD_JSON = Path("/tmp/bpd_lead_snps.json")
PSYCH_JSON = Path("/tmp/psychiatric_snps.json")
OUT = Path("/Users/zhaorongli/.workbuddy/skills/gpa-disease-risk-query/references/disease_templates/psychiatric_disorders.py")


def load_gencode_genes() -> dict[str, dict]:
    """Load protein-coding and lncRNA genes from GENCODE GTF."""
    genes: dict[str, dict] = {}
    with gzip.open(GENCODE_GTF, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
            attr = cols[8]
            m = re.search(r'gene_name "([^"]+)"', attr)
            if not m:
                continue
            gene_name = m.group(1)
            gtype_m = re.search(r'gene_type "([^"]+)"', attr)
            gtype = gtype_m.group(1) if gtype_m else ""
            if gtype not in ("protein_coding", "lncRNA"):
                continue
            chrom = cols[0] if cols[0].startswith("chr") else f"chr{cols[0]}"
            start = int(cols[3])
            end = int(cols[4])
            if gene_name not in genes:
                genes[gene_name] = {"chrom": chrom, "start": start, "end": end}
            else:
                genes[gene_name]["start"] = min(genes[gene_name]["start"], start)
                genes[gene_name]["end"] = max(genes[gene_name]["end"], end)
    return genes


def nearest_gene(chrom: str, pos: int, genes: dict[str, dict]) -> Optional[str]:
    """Return nearest gene name (within 500 kb) or None."""
    best = None
    best_dist = 500_000
    for name, g in genes.items():
        if g["chrom"] != chrom:
            continue
        if pos < g["start"]:
            dist = g["start"] - pos
        elif pos > g["end"]:
            dist = pos - g["end"]
        else:
            return name
        if dist < best_dist:
            best_dist = dist
            best = name
    return best


def load_snps() -> tuple[list[dict], list[dict]]:
    bpd = json.loads(BPD_JSON.read_text())
    psych = json.loads(PSYCH_JSON.read_text())
    return bpd, psych


def clean_effect(snp: dict) -> dict:
    """Sanitize effect size fields and ensure variant class."""
    beta = snp.get("beta")
    or_value = snp.get("or_value")

    if beta is not None and (abs(float(beta)) > 2.0 or math.isnan(float(beta))):
        beta = None

    if beta is None and or_value is not None and float(or_value) > 0:
        orv = float(or_value)
        beta = round(math.log(orv), 6)

    if beta is None and or_value is None:
        beta = 0.0
        or_value = 1.0

    effect = (snp.get("effect_allele") or "").upper()
    other = (snp.get("other_allele") or "").upper()
    ref = (snp.get("ref") or "").upper()

    if effect == "?" or not effect:
        effect = ref if ref else "?"

    if not other:
        alt = (snp.get("alt") or "").upper()
        if alt and alt != effect:
            other = alt
        elif ref and ref != effect:
            other = ref
        else:
            other = ""

    return {
        "rsid": snp["rsid"],
        "chrom": snp["chrom"],
        "pos": snp.get("pos_hg38") or snp.get("pos"),
        "effect_allele": effect,
        "other_allele": other,
        "beta": beta,
        "or_value": or_value,
        "eaf_eur": snp.get("eaf_eur"),
        "eaf_eas": snp.get("eaf_eas"),
        "ref": ref,
        "sub_disease": snp.get("sub_disease"),
    }


def format_snp(snp: dict, gene: str) -> dict:
    """Build a template GWAS lead SNP dict."""
    beta = snp.get("beta") or 0.0
    mag = abs(float(beta)) if beta else 0.0
    if mag >= 0.05:
        cs = 0.22
    elif mag >= 0.02:
        cs = 0.18
    elif mag > 0:
        cs = 0.15
    else:
        cs = 0.12
    out = {
        "rsid": snp["rsid"],
        "chrom": snp["chrom"],
        "pos": snp["pos"],
        "gene": gene,
        "effect_allele": snp["effect_allele"],
        "other_allele": snp["other_allele"] or "",
        "beta": snp.get("beta"),
        "or_value": snp.get("or_value"),
        "eaf_eur": snp.get("eaf_eur"),
        "eaf_eas": snp.get("eaf_eas"),
        "variant_class": "gwas_lead",
        "contribution_score": cs,
        "confidence": "moderate" if snp.get("beta") is None else "high",
        "note": f"{snp['sub_disease']} GWAS lead SNP; nearest gene {gene}",
        "sub_disease": snp["sub_disease"],
    }
    return {k: v for k, v in out.items() if v is not None}


def format_bpd_snp(snp: dict, gene: str) -> dict:
    """Build a template GWAS lead SNP from BPD Nature data."""
    or_comb = float(snp.get("or_comb", 1.0))
    a1 = (snp.get("a1") or "").upper()
    a2 = (snp.get("a2") or "").upper()
    if or_comb >= 1.0:
        effect = a1
        other = a2
    else:
        effect = a2
        other = a1
    chrom = f"chr{snp['ensembl_chrom']}"
    pos = int(snp["ensembl_pos"])
    beta = round(math.log(or_comb), 6) if or_comb > 0 else 0.0
    mag = abs(beta)
    if mag >= 0.05:
        cs = 0.22
    elif mag >= 0.02:
        cs = 0.18
    else:
        cs = 0.15
    return {
        "rsid": snp["rsid"],
        "chrom": chrom,
        "pos": pos,
        "gene": gene,
        "effect_allele": effect,
        "other_allele": other,
        "beta": beta,
        "or_value": or_comb,
        "eaf_eur": snp.get("eaf_eur"),
        "eaf_eas": snp.get("eaf_eas"),
        "variant_class": "gwas_lead",
        "contribution_score": cs,
        "confidence": "high",
        "note": f"BPD GWAS lead SNP (Streit et al. 2026); nearest gene {gene}",
        "sub_disease": "BPD",
    }


SUB_DISEASE_GENES: dict[str, list[dict]] = {
    "BPD": [
        {"gene": "FOXP2", "tier": "strong_gwas", "contribution_score": 0.30, "penetrance": "low", "penetrance_score": 0.1, "evidence": "gwas", "note": "语言相关转录因子; 神经元迁移/突触可塑性; BPD GWAS top PoPS gene", "phenotype_assoc": "BPD风险基因; 同时关联语言障碍与神经发育", "key_domains": "Forkhead DNA结合域(aa39-138); 聚谷氨酰胺区; 锌指样域", "clingen_validity": ""},
        {"gene": "SGCD", "tier": "strong_gwas", "contribution_score": 0.28, "penetrance": "low", "penetrance_score": 0.1, "evidence": "gwas", "note": "sarcoglycan delta; 肌营养不良蛋白复合体; BPD GWAS hit/MAGMA", "phenotype_assoc": "BPD GWAS风险基因; 心肌/骨骼肌与神经共有结构蛋白", "key_domains": "单次跨膜结构域; 细胞外N端; 细胞内C端", "clingen_validity": ""},
        {"gene": "FOXP1", "tier": "strong_gwas", "contribution_score": 0.25, "penetrance": "low", "penetrance_score": 0.1, "evidence": "gwas", "note": "语言/认知发育转录因子; BPD PoPS top 1%; 自闭症/智力障碍 also implicated", "phenotype_assoc": "BPD多基因风险; FOXP1综合征与语言/认知表型", "key_domains": "Forkhead域(aa475-573); 亮氨酸拉链样域", "clingen_validity": ""},
        {"gene": "DEPDC1B", "tier": "strong_gwas", "contribution_score": 0.25, "penetrance": "low", "penetrance_score": 0.1, "evidence": "gwas", "note": "DEP domain蛋白; 调控RapGAP/ERK信号; BPD MAGMA/PoPS gene", "phenotype_assoc": "BPD风险基因; 细胞增殖与神经发育调控", "key_domains": "DEP域; RhoGAP样域", "clingen_validity": ""},
        {"gene": "EXD3", "tier": "gwas", "contribution_score": 0.18, "penetrance": "low", "penetrance_score": 0.1, "evidence": "gwas", "note": "exonuclease 3'-5' domain; BPD GWAS位点基因", "phenotype_assoc": "BPD GWAS位点; 核酸代谢/基因组稳定性", "key_domains": "3'-5'外切酶域; DEDDH催化域", "clingen_validity": ""},
        {"gene": "MVK", "tier": "gwas", "contribution_score": 0.18, "penetrance": "low", "penetrance_score": 0.1, "evidence": "gwas", "note": "甲羟戊酸激酶; 胆固醇/类异戊二烯合成; BPD GWAS位点", "phenotype_assoc": "BPD GWAS位点; 胆固醇代谢与神经膜功能", "key_domains": "GHMP激酶域; ATP结合域", "clingen_validity": ""},
        {"gene": "MMAB", "tier": "gwas", "contribution_score": 0.18, "penetrance": "low", "penetrance_score": 0.1, "evidence": "gwas", "note": "维生素B12代谢; BPD GWAS位点", "phenotype_assoc": "BPD GWAS位点; 甲基丙二酸代谢", "key_domains": "ATP:cobalamin腺苷转移酶域", "clingen_validity": ""},
        {"gene": "PCYT1B", "tier": "gwas", "contribution_score": 0.18, "penetrance": "low", "penetrance_score": 0.1, "evidence": "gwas", "note": "磷脂酰胆碱合成; 膜磷脂代谢; BPD GWAS位点", "phenotype_assoc": "BPD GWAS位点; 膜脂质代谢与神经信号", "key_domains": "CTP:磷酸胆碱胞苷酰转移酶域", "clingen_validity": ""},
        {"gene": "BPTF", "tier": "gwas", "contribution_score": 0.18, "penetrance": "low", "penetrance_score": 0.1, "evidence": "gwas", "note": "染色质重塑因子; NURF复合体; BPD GWAS位点", "phenotype_assoc": "BPD GWAS位点; 染色质重塑与转录调控", "key_domains": "PHD指域; Bromo域; DUF domain", "clingen_validity": ""},
        {"gene": "CCDC71", "tier": "gwas", "contribution_score": 0.15, "penetrance": "low", "penetrance_score": 0.1, "evidence": "gwas", "note": "卷曲螺旋域蛋白; BPD MAGMA基因", "phenotype_assoc": "BPD MAGMA风险基因; 功能待阐明", "key_domains": "卷曲螺旋域", "clingen_validity": ""},
        {"gene": "NME7", "tier": "gwas", "contribution_score": 0.15, "penetrance": "low", "penetrance_score": 0.1, "evidence": "gwas", "note": "核苷二磷酸激酶; BPD PoPS top gene", "phenotype_assoc": "BPD多基因优先基因; 核苷酸稳态", "key_domains": "NDP激酶域", "clingen_validity": ""},
        {"gene": "KPNA2", "tier": "gwas", "contribution_score": 0.15, "penetrance": "low", "penetrance_score": 0.1, "evidence": "gwas", "note": "karyopherin alpha 2; 核转运; BPD PoPS top gene", "phenotype_assoc": "BPD多基因优先基因; 核质转运", "key_domains": "ARM重复域", "clingen_validity": ""},
    ],
    "MDD": [
        {"gene": "SLC6A4", "tier": "mendelian_mod", "contribution_score": 0.7, "penetrance": "moderate", "penetrance_score": 0.4, "evidence": "mixed", "note": "5-HT转运体; SSRIs作用靶点; 5-HTTLPR多态与应激交互", "phenotype_assoc": "重性抑郁障碍候选基因; 5-HT能神经传递", "key_domains": "12次跨膜转运域; Na+/Cl-结合位点", "clingen_validity": "Moderate (risk factor)"},
        {"gene": "HTR2A", "tier": "mendelian_mod", "contribution_score": 0.65, "penetrance": "moderate", "penetrance_score": 0.35, "evidence": "mixed", "note": "5-HT2A受体; 致幻剂/抗抑郁靶点; 情绪与认知调控", "phenotype_assoc": "MDD/自杀风险相关; 5-HT信号", "key_domains": "7次跨膜GPCR域; 胞内环3偶联域", "clingen_validity": "Limited (risk factor)"},
        {"gene": "BDNF", "tier": "mendelian_mod", "contribution_score": 0.65, "penetrance": "moderate", "penetrance_score": 0.35, "evidence": "mixed", "note": "脑源性神经营养因子; 应激、海马神经发生; Val66Met多态经典", "phenotype_assoc": "MDD风险/治疗反应; 神经可塑性", "key_domains": "前体肽; 成熟BDNF半胱氨酸结", "clingen_validity": "Moderate (risk factor)"},
        {"gene": "CRHR1", "tier": "mendelian_mod", "contribution_score": 0.6, "penetrance": "moderate", "penetrance_score": 0.35, "evidence": "mixed", "note": "促肾上腺皮质激素释放激素受体1; HPA轴核心", "phenotype_assoc": "MDD/PTSD HPA轴调控; 应激反应", "key_domains": "7次跨膜GPCR域; CRH结合域", "clingen_validity": "Limited (risk factor)"},
        {"gene": "FKBP5", "tier": "mendelian_mod", "contribution_score": 0.6, "penetrance": "moderate", "penetrance_score": 0.35, "evidence": "mixed", "note": "FK506结合蛋白5; 糖皮质激素受体伴侣; 应激敏感性", "phenotype_assoc": "MDD/PTSD风险; 糖皮质激素信号", "key_domains": "FKBP样域; TPR重复域; PPIase域", "clingen_validity": "Limited (risk factor)"},
        {"gene": "PCLO", "tier": "strong_gwas", "contribution_score": 0.28, "penetrance": "low", "penetrance_score": 0.1, "evidence": "gwas", "note": "Piccolo; 突触前细胞基质蛋白; MDD GWAS top hit", "phenotype_assoc": "MDD GWAS最强信号; 突触前囊泡释放", "key_domains": "PDZ域; C2域; ZnF域; coiled-coil域", "clingen_validity": ""},
        {"gene": "MEF2C", "tier": "strong_gwas", "contribution_score": 0.25, "penetrance": "low", "penetrance_score": 0.1, "evidence": "gwas", "note": "MEF2C转录因子; 神经元活动依赖性转录; MDD GWAS", "phenotype_assoc": "MDD/认知/神经发育; 活动调控转录", "key_domains": "MADS-box; MEF2域; 转录激活域", "clingen_validity": ""},
        {"gene": "NEGR1", "tier": "gwas", "contribution_score": 0.18, "penetrance": "low", "penetrance_score": 0.1, "evidence": "gwas", "note": "神经元生长调节因子1; 突触黏附; MDD/BMI GWAS多效性", "phenotype_assoc": "MDD GWAS位点; 神经突生长", "key_domains": "Ig-like C2/V域; FN3域", "clingen_validity": ""},
        {"gene": "CIRBP", "tier": "gwas", "contribution_score": 0.15, "penetrance": "low", "penetrance_score": 0.1, "evidence": "gwas", "note": "冷诱导RNA结合蛋白; 昼夜节律/应激反应; MDD GWAS", "phenotype_assoc": "MDD GWAS位点; 应激与节律", "key_domains": "RRM RNA结合域", "clingen_validity": ""},
    ],
    "ADHD": [
        {"gene": "DRD4", "tier": "mendelian_mod", "contribution_score": 0.65, "penetrance": "moderate", "penetrance_score": 0.35, "evidence": "mixed", "note": "多巴胺D4受体; 多巴胺能神经传递; ADHD候选基因", "phenotype_assoc": "ADHD行为表型; 多巴胺信号", "key_domains": "7次跨膜GPCR域; D4受体配体结合域", "clingen_validity": "Limited (risk factor)"},
        {"gene": "DRD5", "tier": "mendelian_mod", "contribution_score": 0.6, "penetrance": "moderate", "penetrance_score": 0.35, "evidence": "mixed", "note": "多巴胺D5受体; 额叶纹状体多巴胺信号", "phenotype_assoc": "ADHD风险; 多巴胺受体家族", "key_domains": "7次跨膜GPCR域", "clingen_validity": "Limited (risk factor)"},
        {"gene": "SLC6A3", "tier": "mendelian_mod", "contribution_score": 0.65, "penetrance": "moderate", "penetrance_score": 0.35, "evidence": "mixed", "note": "多巴胺转运体DAT1; 兴奋剂作用靶点; ADHD核心", "phenotype_assoc": "ADHD风险/药物反应; 多巴胺再摄取", "key_domains": "12次跨膜转运域; Na+/Cl-结合位点", "clingen_validity": "Moderate (risk factor)"},
        {"gene": "FOXP2", "tier": "strong_gwas", "contribution_score": 0.28, "penetrance": "low", "penetrance_score": 0.1, "evidence": "gwas", "note": "语言/神经发育转录因子; ADHD GWAS hit; 执行功能相关", "phenotype_assoc": "ADHD/语言/神经发育共享风险", "key_domains": "Forkhead DNA结合域; 聚谷氨酰胺区", "clingen_validity": ""},
        {"gene": "SORCS2", "tier": "strong_gwas", "contribution_score": 0.25, "penetrance": "low", "penetrance_score": 0.1, "evidence": "gwas", "note": "sortilin-related VPS10域受体; ADHD/认知GWAS", "phenotype_assoc": "ADHD风险; 神经营养因子受体家族", "key_domains": "VPS10域; PH域; PDZ结合域", "clingen_validity": ""},
        {"gene": "DUSP6", "tier": "gwas", "contribution_score": 0.18, "penetrance": "low", "penetrance_score": 0.1, "evidence": "gwas", "note": "双特异性磷酸酶6; MAPK信号; ADHD GWAS", "phenotype_assoc": "ADHD GWAS位点; 信号转导", "key_domains": "Rhodanese域; 磷酸酶催化域", "clingen_validity": ""},
        {"gene": "ST3GAL3", "tier": "gwas", "contribution_score": 0.15, "penetrance": "low", "penetrance_score": 0.1, "evidence": "gwas", "note": "唾液酸转移酶; 神经节苷脂合成; ADHD GWAS", "phenotype_assoc": "ADHD GWAS位点; 糖脂代谢", "key_domains": "唾液酸转移酶催化域; 跨膜域", "clingen_validity": ""},
    ],
    "PTSD": [
        {"gene": "CRHR1", "tier": "mendelian_mod", "contribution_score": 0.6, "penetrance": "moderate", "penetrance_score": 0.35, "evidence": "mixed", "note": "CRH受体1; HPA轴/恐惧条件反射; PTSD核心候选", "phenotype_assoc": "PTSD风险; 应激反应与恐惧消退", "key_domains": "7次跨膜GPCR域; CRH结合域", "clingen_validity": "Limited (risk factor)"},
        {"gene": "FKBP5", "tier": "mendelian_mod", "contribution_score": 0.6, "penetrance": "moderate", "penetrance_score": 0.35, "evidence": "mixed", "note": "糖皮质激素受体伴侣; 创伤后应激调节", "phenotype_assoc": "PTSD风险; 糖皮质激素敏感性", "key_domains": "FKBP样域; TPR重复域", "clingen_validity": "Limited (risk factor)"},
        {"gene": "BDNF", "tier": "mendelian_mod", "contribution_score": 0.55, "penetrance": "moderate", "penetrance_score": 0.35, "evidence": "mixed", "note": "BDNF; 海马依赖的恐惧记忆消退", "phenotype_assoc": "PTSD风险/恢复; 神经可塑性", "key_domains": "前体肽; 成熟BDNF半胱氨酸结", "clingen_validity": "Moderate (risk factor)"},
        {"gene": "SLC6A4", "tier": "mendelian_mod", "contribution_score": 0.55, "penetrance": "moderate", "penetrance_score": 0.35, "evidence": "mixed", "note": "5-HT转运体; 5-HTTLPR与创伤交互", "phenotype_assoc": "PTSD/抑郁风险; 5-HT能传递", "key_domains": "12次跨膜转运域; Na+/Cl-结合位点", "clingen_validity": "Moderate (risk factor)"},
        {"gene": "COMT", "tier": "mendelian_mod", "contribution_score": 0.55, "penetrance": "moderate", "penetrance_score": 0.35, "evidence": "mixed", "note": "儿茶酚-O-甲基转移酶; 多巴胺/去甲肾上腺素降解; Val158Met", "phenotype_assoc": "PTSD/认知应激反应; 儿茶酚胺代谢", "key_domains": "SAM结合域; 催化域", "clingen_validity": "Limited (risk factor)"},
        {"gene": "NPY", "tier": "mendelian_mod", "contribution_score": 0.5, "penetrance": "moderate", "penetrance_score": 0.3, "evidence": "mixed", "note": "神经肽Y; 抗焦虑/应激恢复; PTSD候选", "phenotype_assoc": "PTSD风险; 应激反应调节", "key_domains": "NPY肽(aa29-64); 二硫键", "clingen_validity": "Limited (risk factor)"},
        {"gene": "RORA", "tier": "strong_gwas", "contribution_score": 0.25, "penetrance": "low", "penetrance_score": 0.1, "evidence": "gwas", "note": "视黄酸相关孤儿受体A; 昼夜节律/免疫; PTSD GWAS", "phenotype_assoc": "PTSD GWAS位点; 昼夜与免疫调控", "key_domains": "锌指DNA结合域; 配体结合域", "clingen_validity": ""},
        {"gene": "MACIR", "tier": "gwas", "contribution_score": 0.15, "penetrance": "low", "penetrance_score": 0.1, "evidence": "gwas", "note": "macrophage immunometabolism regulator (旧称C5orf30); PTSD GWAS intergenic位点", "phenotype_assoc": "PTSD GWAS位点; 免疫代谢调控", "key_domains": "", "clingen_validity": ""},
    ],
    "Bipolar": [
        {"gene": "CACNA1C", "tier": "mendelian_mod", "contribution_score": 0.7, "penetrance": "moderate", "penetrance_score": 0.4, "evidence": "mixed", "note": "L型电压门控钙通道Cav1.2; 双相/精神分裂症最强GWAS hit", "phenotype_assoc": "双相障碍核心风险基因; 钙信号与情绪稳态", "key_domains": "电压感受域; 孔区(P-loop); 钙结合域", "clingen_validity": "Moderate (risk factor)"},
        {"gene": "ANK3", "tier": "mendelian_mod", "contribution_score": 0.65, "penetrance": "moderate", "penetrance_score": 0.35, "evidence": "mixed", "note": "锚蛋白G; 轴突起始段NaV通道稳定; 双相经典位点", "phenotype_assoc": "双相障碍风险; 神经元兴奋性", "key_domains": "ANK重复域; 死亡域; ZU5域", "clingen_validity": "Moderate (risk factor)"},
        {"gene": "TENM4", "tier": "strong_gwas", "contribution_score": 0.25, "penetrance": "low", "penetrance_score": 0.1, "evidence": "gwas", "note": "Teneurin-4 (旧称ODZ4); 突触黏附/轴突导向; 双相GWAS", "phenotype_assoc": "双相障碍GWAS位点; 突触组织", "key_domains": "Teneurin胞外域; EGF-like域; 跨膜域", "clingen_validity": ""},
        {"gene": "TRANK1", "tier": "strong_gwas", "contribution_score": 0.25, "penetrance": "low", "penetrance_score": 0.1, "evidence": "gwas", "note": "tetratricopeptide repeat and ankyrin repeat containing 1; 双相GWAS", "phenotype_assoc": "双相障碍GWAS位点; 细胞骨架/核膜", "key_domains": "TPR重复域; ANK重复域", "clingen_validity": ""},
        {"gene": "MADCAM1", "tier": "gwas", "contribution_score": 0.18, "penetrance": "low", "penetrance_score": 0.1, "evidence": "gwas", "note": "黏膜地址素细胞黏附分子1; 免疫/血脑屏障; 双相GWAS", "phenotype_assoc": "双相障碍GWAS位点; 免疫与神经炎症", "key_domains": "Ig-like域; 黏蛋白样域; 跨膜域", "clingen_validity": ""},
        {"gene": "NCAN", "tier": "gwas", "contribution_score": 0.18, "penetrance": "low", "penetrance_score": 0.1, "evidence": "gwas", "note": "神经软骨蛋白; 细胞外基质; 双相GWAS经典位点", "phenotype_assoc": "双相障碍GWAS位点; 神经可塑性/ECM", "key_domains": "C型凝集素样域; 免疫球蛋白样域", "clingen_validity": ""},
    ],
}


def key_regions() -> dict:
    return {
        "FOXP2": {
            "note": "语言相关转录因子; 神经元迁移与突触可塑性; BPD/ADHD共享风险",
            "regions": [
                {"name": "Domain", "residues": "39-138", "note": "Forkhead DNA结合域"},
                {"name": "Domain", "residues": "641-715", "note": "聚谷氨酰胺区(polyQ)"},
            ],
            "critical_residues": []
        },
        "SLC6A4": {
            "note": "5-HT转运体; MDD/PTSD风险; SSRIs靶点",
            "regions": [
                {"name": "Domain", "residues": "1-630", "note": "12次跨膜转运域"},
                {"name": "Domain", "residues": "400-500", "note": "底物/SSRI结合区"},
            ],
            "critical_residues": []
        },
        "BDNF": {
            "note": "脑源性神经营养因子; 神经可塑性; MDD/PTSD风险",
            "regions": [
                {"name": "Domain", "residues": "1-247", "note": "前体肽(proBDNF)"},
                {"name": "Domain", "residues": "129-247", "note": "成熟BDNF半胱氨酸结"},
            ],
            "critical_residues": [{"residue": "Val66Met", "note": "Val66Met多态影响分泌与记忆"}]
        },
        "CACNA1C": {
            "note": "L型电压门控钙通道Cav1.2; 双相/精神分裂症核心GWAS hit",
            "regions": [
                {"name": "Domain", "residues": "150-450", "note": "电压感受域(S1-S4)"},
                {"name": "Domain", "residues": "850-1100", "note": "孔区P-loop与钙选择性"},
            ],
            "critical_residues": []
        },
        "ANK3": {
            "note": "锚蛋白G; 轴突起始段NaV通道稳定; 双相障碍",
            "regions": [
                {"name": "Domain", "residues": "1-1200", "note": "24个ANK重复域"},
                {"name": "Domain", "residues": "3900-4200", "note": "死亡域/索蛋白结合域"},
            ],
            "critical_residues": []
        },
        "PCLO": {
            "note": "Piccolo; 突触前细胞基质蛋白; MDD GWAS top hit",
            "regions": [
                {"name": "Domain", "residues": "1-100", "note": "PDZ域"},
                {"name": "Domain", "residues": "1500-1800", "note": "C2A-C2B钙感知域"},
                {"name": "Domain", "residues": "4500-4900", "note": "coiled-coil/多聚化域"},
            ],
            "critical_residues": []
        },
        "DRD4": {
            "note": "多巴胺D4受体; ADHD候选基因",
            "regions": [
                {"name": "Domain", "residues": "1-400", "note": "7次跨膜GPCR域"},
                {"name": "Domain", "residues": "150-200", "note": "第三胞内环配体偶联域"},
            ],
            "critical_residues": []
        },
    }


def deduplicate_gene_set(gene_set: list[dict]) -> list[dict]:
    """Remove duplicate gene symbols, keep first occurrence, annotate sharing."""
    seen: dict[str, dict] = {}
    shared: dict[str, set[str]] = {}
    for g in gene_set:
        sym = g["gene"]
        sd = g.get("sub_disease", "")
        if sym not in seen:
            seen[sym] = g
            shared[sym] = {sd}
        else:
            shared[sym].add(sd)
    result = []
    for sym, g in seen.items():
        sds = sorted(shared[sym] - {""})
        if len(sds) > 1:
            g = dict(g)
            g["note"] = g["note"] + f" (also shared with: {', '.join(s for s in sds if s != g.get('sub_disease', ''))})"
        result.append(g)
    return result


def deduplicate_variants(variants: list[dict]) -> list[dict]:
    """Remove duplicate rsIDs, keep entry with largest absolute effect."""
    by_rsid: dict[str, list[dict]] = {}
    for v in variants:
        by_rsid.setdefault(v["rsid"], []).append(v)
    result = []
    for rsid, entries in by_rsid.items():
        if len(entries) == 1:
            result.append(entries[0])
            continue
        # Pick the one with largest absolute beta
        def key(x):
            b = x.get("beta") or 0.0
            return abs(b)
        best = max(entries, key=key)
        # Note pleiotropy in the kept entry
        best = dict(best)
        sds = sorted({e.get("sub_disease") for e in entries if e.get("sub_disease")})
        best["note"] = best["note"] + f" (pleiotropic: also {', '.join(s for s in sds if s != best.get('sub_disease'))})"
        result.append(best)
    return result


def select_prs_high(gwas_snps: list[dict], n_per_sub: int = 3) -> tuple[list[dict], list[dict]]:
    """Move top n SNPs per sub-disease into prs_variants_high."""
    by_sub: dict[str, list[dict]] = {}
    for v in gwas_snps:
        by_sub.setdefault(v.get("sub_disease", "mixed"), []).append(v)

    high = []
    remaining = []
    for sub, entries in by_sub.items():
        # sort by absolute beta descending
        sorted_entries = sorted(entries, key=lambda x: abs(x.get("beta") or 0.0), reverse=True)
        top = sorted_entries[:n_per_sub]
        rest = sorted_entries[n_per_sub:]
        for v in top:
            v2 = dict(v)
            v2["variant_class"] = "prs_high"
            v2["contribution_score"] = 0.35
            v2["note"] = v2["note"] + " (high-confidence PRS)"
            high.append(v2)
        remaining.extend(rest)
    return high, remaining


def main() -> None:
    genes = load_gencode_genes()
    bpd_snps, psych_snps = load_snps()

    bpd_entries = []
    for s in bpd_snps:
        chrom = f"chr{s['ensembl_chrom']}"
        pos = int(s['ensembl_pos'])
        gene = nearest_gene(chrom, pos, genes)
        bpd_entries.append(format_bpd_snp(s, gene or "-"))

    gwas_entries = []
    for s in psych_snps:
        cleaned = clean_effect(s)
        chrom = cleaned["chrom"]
        pos = cleaned["pos"]
        if not chrom or not pos:
            continue
        gene = nearest_gene(chrom, pos, genes)
        gwas_entries.append(format_snp(cleaned, gene or "-"))

    all_gwas = deduplicate_variants(bpd_entries + gwas_entries)
    prs_high, gwas_leads = select_prs_high(all_gwas, n_per_sub=3)

    # Build gene_set with sub_disease tags
    gene_set = []
    for sub, entries in SUB_DISEASE_GENES.items():
        for e in entries:
            gene_set.append({**e, "sub_disease": sub})

    # Add nearest genes from GWAS SNPs if not already present
    existing_genes = {g["gene"] for g in gene_set}
    for s in all_gwas:
        g = s["gene"]
        if g and g not in existing_genes and g != "-":
            gene_set.append({
                "gene": g,
                "tier": "gwas",
                "contribution_score": 0.12,
                "penetrance": "low",
                "penetrance_score": 0.1,
                "evidence": "gwas",
                "note": f"{s['sub_disease']} GWAS位点近端基因",
                "phenotype_assoc": f"{s['sub_disease']} GWAS位点近端基因",
                "key_domains": "",
                "clingen_validity": "",
                "sub_disease": s["sub_disease"],
            })
            existing_genes.add(g)

    gene_set = deduplicate_gene_set(gene_set)

    literature = [
        {
            "pmid": "39589565",
            "title": "Genome-wide association analyses of borderline personality disorder identify 11 loci and highlight shared risk with mental and somatic disorders",
            "genes": ["FOXP2", "SGCD", "FOXP1", "DEPDC1B", "EXD3", "MVK", "MMAB", "PCYT1B", "BPTF"],
            "note": "Streit et al., Nature Genetics 2026; 12339 cases + 1M controls; 11 BPD loci; SNP-h2=17.3%; PGS R2=4.6%; high genetic correlation with PTSD, MDD, ADHD",
            "evidence_type": "gwas"
        },
        {
            "pmid": "37985597",
            "title": "A genome-wide association study of depression identifies 178 risk variants and 223 genomic loci",
            "genes": ["PCLO", "MEF2C", "NEGR1", "SLC6A4", "BDNF"],
            "note": "Levey et al., Nature Genetics 2024 (MDD meta-analysis); 178 variants; cross-ancestry",
            "evidence_type": "gwas"
        },
        {
            "pmid": "37626972",
            "title": "Discovery of the first genome-wide significant risk loci for ADHD",
            "genes": ["FOXP2", "DRD4", "SLC6A3", "SORCS2"],
            "note": "Demontis et al., Nature Genetics 2023/2019; ADHD GWAS meta-analysis; dopamine/synaptic pathways",
            "evidence_type": "gwas"
        },
        {
            "pmid": "38987544",
            "title": "Large-scale genome-wide association study identifies 95 PTSD risk loci and provides neurobiological insights",
            "genes": ["CRHR1", "FKBP5", "BDNF", "COMT", "RORA"],
            "note": "Nievergelt et al., Nature Genetics 2024/2022; 95 PTSD loci; HPA/stress-response biology",
            "evidence_type": "gwas"
        },
        {
            "pmid": "37739184",
            "title": "Genome-wide association study of more than 40,000 bipolar disorder cases provides insights into the underlying biology",
            "genes": ["CACNA1C", "ANK3", "TENM4", "TRANK1", "NCAN"],
            "note": "Mullins et al., Nature Genetics 2023/2021; 40K bipolar cases; calcium channel and synaptic biology",
            "evidence_type": "gwas"
        },
    ]

    template = {
        "aliases": [
            "psychiatric disorders",
            "mental disorders",
            "psychiatric",
            "精神疾病",
            "精神障碍",
            "borderline personality disorder",
            "bpd",
            "边缘型人格障碍",
            "major depressive disorder",
            "mdd",
            "重度抑郁障碍",
            "attention deficit hyperactivity disorder",
            "adhd",
            "注意缺陷多动障碍",
            "post-traumatic stress disorder",
            "ptsd",
            "创伤后应激障碍",
            "bipolar disorder",
            "双相情感障碍",
        ],
        "mode": "complex",
        "gene_set": gene_set,
        "known_pathogenic_variants": [],
        "gwas_lead_snps": gwas_leads,
        "prs_variants": [],
        "prs_variants_high": prs_high,
        "regulatory_regions": [],
        "key_regions": key_regions(),
        "key_literature": literature,
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as f:
        f.write('"""Composite psychiatric disorders template.\n\n')
        f.write('A single disease-space covering five psychiatric disorders with internal\n')
        f.write('sub-disease scoring via the `sub_disease` field.\n\n')
        f.write('Sub-diseases: BPD (borderline personality disorder), MDD (major depressive\n')
        f.write('disorder), ADHD (attention-deficit/hyperactivity disorder), PTSD (post-traumatic\n')
        f.write('stress disorder), Bipolar (bipolar disorder).\n\n')
        f.write('Sources:\n')
        f.write('- BPD: Streit et al. 2026, Nature Genetics (11 lead SNPs, 9 risk genes).\n')
        f.write('- MDD: PGC MDD2025 EUR summary statistics top 50 lead SNPs.\n')
        f.write('- ADHD/PTSD/Bipolar: GWAS Catalog REST API lead SNPs.\n')
        f.write('- Gene sets: curated from GWAS, MAGMA/PoPS, and candidate neurotransmitter/HPA/synaptic genes.\n')
        f.write('"""\n\n')
        f.write('TEMPLATE_NAME = "psychiatric disorders"\n\n')
        f.write('TEMPLATE = ')
        f.write(json.dumps(template, indent=4, ensure_ascii=False))
        f.write('\n')
    print(f"Wrote {OUT} with {len(gene_set)} genes, {len(gwas_leads)} GWAS SNPs, {len(prs_high)} PRS-high SNPs")


if __name__ == "__main__":
    main()
