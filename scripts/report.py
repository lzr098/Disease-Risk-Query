"""Markdown report generation for disease risk query."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Any, Optional

from constants import COMPLEX_WEIGHT_GWAS_COMMON


def _get(v: dict, *keys: str) -> Any:
    """Get first matching key from variant dict (handles GPA lowercase output)."""
    for k in keys:
        if k in v:
            return v[k]
    return ""


def _variant_table_row(v: dict) -> str:
    gene = _get(v, "GENE", "gene")
    chrom = _get(v, "CHROM", "chrom")
    pos = _get(v, "POS", "pos")
    ref = _get(v, "REF", "ref")
    alt = _get(v, "ALT", "alt")
    hgvsp = _get(v, "HGVSp", "hgvsp", "primary_hgvsp")
    hgvsc = _get(v, "HGVSc", "hgvsc", "primary_hgvsc")
    clin_sig = _get(v, "CLIN_SIG", "clinvar")
    impact = _get(v, "IMPACT", "impact", "primary_impact")
    cons = _get(v, "Consequence", "consequence", "primary_consequence")
    gt = _get(v, "GT", "gt")
    gnomad = _get(v, "gnomAD_AF", "gnomad_af")
    flag_parts = list(v.get("_drq_flags", []))
    if v.get("gwas_proximal"):
        flag_parts.append(f"GWAS:{v.get('nearest_gwas_snp', '')}")
    if v.get("literature_support"):
        flag_parts.append("LIT")
    flags = ", ".join(flag_parts)

    # Template-based contribution/penetrance annotations
    gene_contrib = v.get("_gene_contribution")
    gene_pen = v.get("_gene_penetrance")
    contrib_str = f"{gene_contrib}" if gene_contrib is not None else "-"
    pen_str = f"{gene_pen}" if gene_pen else "-"

    return (
        f"| {gene} | `{chrom}:{pos}:{ref}>{alt}` | {hgvsc} | {hgvsp} | "
        f"{clin_sig} | {impact}/{cons} | {gt} | {gnomad} | {contrib_str} | {pen_str} | {flags} |"
    )


def generate_report(
    disease_name: str,
    hpo_id: Optional[str],
    hpo_name: Optional[str],
    sex: str,
    age: Optional[int],
    gene_set_result: dict,
    gpa_result: dict,
    score_result: dict,
    output_path: Path,
    apoe_result: Optional[dict] = None,
    gwas_summary: Optional[dict] = None,
    literature_summary: Optional[dict] = None,
    disease_reference: Optional[dict] = None,
    gwas_lead_snps: Optional[list[dict]] = None,
    vcf_qc: Optional[dict] = None,
    disease_mode: str = "mendelian",
) -> Path:
    """Generate the final Markdown risk/contribution report."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    is_complex = disease_mode == "complex"

    lines: list[str] = []
    title = f"疾病遗传贡献度评估报告：{disease_name}" if is_complex else f"疾病风险查询报告：{disease_name}"
    lines.append(f"# {title}")
    lines.append("")
    lines.append("## 1. 查询摘要")
    lines.append(f"- **目标疾病**：{disease_name}")
    lines.append(f"- **评估模式**：{'复杂表型遗传贡献度' if is_complex else '孟德尔/致病性风险'}")
    lines.append(f"- **HPO ID**：{hpo_id or 'N/A'} ({hpo_name or 'N/A'})")
    lines.append(f"- **样本性别**：{sex} | **年龄**：{age if age is not None else '未知'}")
    lines.append(f"- **查询日期**：{datetime.now().isoformat(timespec='minutes')}")
    lines.append("")
    lines.append("### VCF 质量与完整性")
    if vcf_qc and vcf_qc.get("checked"):
        lines.append(
            f"- **锚定位点检出率**：{vcf_qc.get('anchor_snps_present', 0)} / "
            f"{vcf_qc.get('anchor_snps_checked', 0)} "
            f"({vcf_qc.get('presence_rate', 0):.0%})"
        )
        lines.append(f"- **VCF 总变异数（近似）**：{vcf_qc.get('total_variants', 'N/A'):,}")
        if vcf_qc.get("common_variants_filtered"):
            lines.append(
                "- **⚠️ 常见变异过滤提示**：该 VCF 的锚定位点检出率低于阈值，"
                "说明常见 SNP 位点未在 VCF 中保留。对于已基因分型的 callset，"
                "未保留的常见 SNP 可推断为 ref/ref（0/0），并非测序漏检；"
                "复杂表型模式下 GWAS 维度已按此情况自动折算。"
            )
        else:
            lines.append("- **VCF 完整性**：锚定位点检出率正常，GWAS 维度可正常评估。")
    else:
        lines.append(f"- **VCF QC**：{vcf_qc.get('note', '未执行') if vcf_qc else '未执行'}")
    lines.append("")

    if apoe_result:
        lines.append("## 2. APOE 基因型")
        if apoe_result.get("present"):
            lines.append(f"- **rs7412**：{apoe_result.get('rs7412_status')}")
            lines.append(f"- **rs429358**：{apoe_result.get('rs429358_status')}")
            inferred = apoe_result.get("inferred_allele") or "无法推断"
            lines.append(f"- **推断 APOE 等位基因**：{inferred}")
        else:
            lines.append("- **APOE 位点未检出**：VCF 在 chr19:44,908,637-44,912,685 区域无变异记录。")
        if apoe_result.get("warning"):
            lines.append(f"- **提示**：{apoe_result['warning']}")
        lines.append("")

    lines.append("## 3. 关联基因集")
    sources = gene_set_result.get("sources", {})
    lines.append(f"- **数据库来源**：OMIM {sources.get('omim', 0)} 个，HPO {sources.get('hpo', 0)} 个")
    lines.append(f"- **文献扩展**：{sources.get('extra', 0)} 个新基因")
    lines.append(f"- **合并去重后基因总数**：{gene_set_result.get('total', 0)}")
    lines.append("")
    lines.append("**基因列表（前 30）**：")
    genes = gene_set_result.get("merged_genes", [])
    lines.append(", ".join(genes[:30]) + (" ..." if len(genes) > 30 else ""))
    lines.append("")

    tier1 = gpa_result.get("tier1_variants", [])
    tier2 = gpa_result.get("tier2_variants", [])
    tier3 = gpa_result.get("tier3_variants", [])
    multi_hit = gpa_result.get("multi_hit", [])

    # Enrich variants with template-based gene contribution / penetrance
    gene_contribution_map = (disease_reference or {}).get("gene_contribution_map", {})
    gene_penetrance_map = (disease_reference or {}).get("gene_penetrance_map", {})
    for v in tier1 + tier2 + tier3:
        gene = _get(v, "GENE", "gene")
        if gene and gene in gene_contribution_map:
            v["_gene_contribution"] = gene_contribution_map[gene]
            v["_gene_penetrance"] = gene_penetrance_map.get(gene, "")

    lines.append("## 4. 变异匹配结果")
    lines.append("")
    if is_complex:
        lines.append(f"- **Tier 1（高外显率）**：{len(tier1)} 个")
        lines.append(f"- **Tier 2（可能相关）**：{len(tier2)} 个")
        lines.append(f"- **Tier 3（罕见功能 / 去噪后）**：{len(tier3)} 个")
        lines.append("")
        lines.append("> 在复杂表型模式下，Tier 仅表示变异置信度层级；最终解释以「遗传贡献度」为准。")
        lines.append("")
    else:
        lines.append(f"### Tier 1（高置信度）— {len(tier1)} 个")
        if tier1:
            lines.append("| 基因 | 位点 | cDNA | 蛋白 | ClinVar | 影响/类型 | 合子性 | gnomAD AF | 基因贡献 | 外显率 | 标志 |")
            lines.append("|------|------|------|------|---------|-----------|--------|-----------|----------|--------|------|")
            for v in tier1:
                lines.append(_variant_table_row(v))
        else:
            lines.append("_未检出_")
        lines.append("")

        lines.append(f"### Tier 2（可能相关）— {len(tier2)} 个")
        if tier2:
            lines.append("| 基因 | 位点 | cDNA | 蛋白 | ClinVar | 影响/类型 | 合子性 | gnomAD AF | 基因贡献 | 外显率 | 标志 |")
            lines.append("|------|------|------|------|---------|-----------|--------|-----------|----------|--------|------|")
            for v in tier2:
                lines.append(_variant_table_row(v))
        else:
            lines.append("_未检出_")
        lines.append("")

    # Tier 3 table shown for both modes (often the only findings in complex traits)
    tier_label = "Tier 3（弱证据 / 去噪后）" if not is_complex else "Tier 3（罕见功能变异 / 去噪后）"
    lines.append(f"### {tier_label} — {len(tier3)} 个")
    if tier3:
        lines.append("| 基因 | 位点 | cDNA | 蛋白 | ClinVar | 影响/类型 | 合子性 | gnomAD AF | 基因贡献 | 外显率 | 标志 |")
        lines.append("|------|------|------|------|---------|-----------|--------|-----------|----------|--------|------|")
        for v in tier3:
            lines.append(_variant_table_row(v))
    else:
        lines.append("_未检出_")
    lines.append("")

    if not is_complex:
        lines.append("### Multi-hit 检测")
        if multi_hit:
            for item in multi_hit:
                if isinstance(item, dict):
                    gene = item.get("gene", "") or "(unknown gene)"
                    count = item.get("variant_count", item.get("count", "?"))
                    phase_result = item.get("phase_result", {})
                    phase = phase_result.get("status", item.get("phase", "unknown"))
                    lines.append(f"- **{gene}**：命中 {count} 个变异，相位状态：{phase}")
                else:
                    lines.append(f"- {item}")
        else:
            lines.append("_未检出同基因多个命中变异_")
        lines.append("")

    lines.append("## 5. GWAS 与文献证据")
    if disease_reference:
        meta = disease_reference.get("metadata", {})
        lines.append(f"- **疾病参考缓存**：{meta.get('source', 'unknown')} / {meta.get('created', 'N/A')}")
        counts = meta.get("counts", {})
        lines.append(f"- **核心基因**：{counts.get('core_genes', 0)} 个")
        lines.append(f"- **已知 ClinVar 致病突变**：{counts.get('clinvar_variants', 0)} 个")
        lines.append(f"- **GWAS lead SNPs**：{counts.get('gwas_snps', 0)} 个")
        lines.append(f"- **关键文献条目**：{counts.get('literature_entries', 0)} 条")
    lines.append("")

    if gwas_lead_snps:
        lines.append("### GWAS lead SNP 直接检出")
        hits = [s for s in gwas_lead_snps if s.get("sample_gt")]
        if hits:
            lines.append(f"- 在样本中直接检出的已知 GWAS lead SNP：**{len(hits)} / {len(gwas_lead_snps)}**")
            lines.append("| SNP | 基因 | 位点 | 样本基因型 | 效应等位 | beta | OR | 贡献 | 说明 |")
            lines.append("|-----|------|------|------------|----------|------|----|------|------|")
            for s in hits:
                gt = s["sample_gt"]
                beta = s.get("beta")
                or_val = s.get("or")
                contrib = s.get("contribution_score")
                beta_str = f"{beta:.4f}" if beta is not None else "-"
                or_str = f"{or_val:.2f}" if or_val is not None else "-"
                contrib_str = f"{contrib:.3f}" if contrib is not None else "-"
                lines.append(
                    f"| {s.get('rsid', '')} | {s.get('gene', '')} | "
                    f"{gt.get('chrom', '')}:{gt.get('pos', '')} | {gt.get('gt', '')} | "
                    f"{s.get('effect_allele', '')} | {beta_str} | {or_str} | {contrib_str} | {s.get('note', '')} |"
                )
        else:
            lines.append(f"- 已知 GWAS lead SNP 均未在 VCF 中保留（共 {len(gwas_lead_snps)} 个）。")
            lines.append("- 在已基因分型的 callset 中，未保留的常见 SNP 应推断为 ref/ref（0/0），而非漏检。")
        lines.append("")

    if gwas_summary:
        lines.append("### GWAS 位点覆盖（±500 kb）")
        lines.append(f"- 本样本中落在 GWAS lead SNP ±500 kb 窗口内的变异：**{gwas_summary.get('hit_count', 0)} 个**")
        hit_genes = gwas_summary.get("hit_genes", [])
        if hit_genes:
            lines.append(f"- 涉及 GWAS 基因：{', '.join(hit_genes)}")
        else:
            lines.append("- 未涉及已知 GWAS 基因座")
        lines.append("")

    if literature_summary:
        lines.append("### 关键文献支持")
        lines.append(f"- 有文献支持的变异：**{len(literature_summary.get('variant_hits', []))} 个**")
        lit_genes = literature_summary.get("gene_hits", [])
        if lit_genes:
            lines.append(f"- 涉及文献基因：{', '.join(lit_genes)}")
        else:
            lines.append("- 未涉及文献支持基因")
        # List top PMIDs
        seen_pmids = set()
        for v in literature_summary.get("variant_hits", []):
            for entry in v.get("literature_support", []):
                pmid = entry.get("pmid")
                if pmid and pmid not in seen_pmids:
                    seen_pmids.add(pmid)
                    title = entry.get("title", "")
                    note = entry.get("note", "")
                    lines.append(f"  - PMID:{pmid} — {note}")
        lines.append("")

    comp = score_result.get("components", {})
    if is_complex:
        lines.append("## 6. 遗传贡献度评估")
        total = score_result.get("total_score", 0)
        level = score_result.get("contribution_level", "无明确遗传贡献")
        meaning = score_result.get("contribution_meaning", "")
        lines.append(f"- **遗传贡献度总分**：{total}/100")
        lines.append(f"- **贡献度等级**：**{level}**")
        lines.append(f"- **含义**：{meaning}")
        lines.append("")
        lines.append("**分数组成**：")
        lines.append(f"- 高外显率单基因变异：{comp.get('monogenic_score', 0)} / 30")
        lines.append(f"- 罕见功能变异：{comp.get('rare_functional_score', 0)} / 25")
        gwas_adj = score_result.get("gwas_adjustment")
        if gwas_adj:
            lines.append(
                f"- GWAS 风险等位基因：{comp.get('gwas_score', 0)} / "
                f"{gwas_adj.get('effective_gwas_weight', COMPLEX_WEIGHT_GWAS_COMMON)} "
                f"（原 {comp.get('gwas_original_score', comp.get('gwas_score', 0))}，"
                f"因 VCF 过滤按 {gwas_adj.get('downweight_factor', 1):.0%} 折算）"
            )
            lines.append(f"- GWAS 未评估权重：{gwas_adj.get('unassessed_gwas_weight', 0)} / {COMPLEX_WEIGHT_GWAS_COMMON}")
            lines.append(f"> {gwas_adj.get('note', '')}")
        else:
            lines.append(f"- GWAS 风险等位基因：{comp.get('gwas_score', 0)} / 25")
        lines.append(f"- 文献/通路支持：{comp.get('literature_score', 0)} / 10")
        lines.append(f"- 人群罕见度：{round(comp.get('rarity_score', 0), 1)} / 10")
        if comp.get("min_gnomad_af") is not None:
            lines.append(f"- 最低 gnomAD AF：{comp['min_gnomad_af']}")
        lines.append("")

        # Per-variant contribution table
        rare_details = comp.get("rare_functional_details", [])
        gwas_hits = comp.get("gwas_hits", [])
        if rare_details or gwas_hits:
            lines.append("### 关键突变/位点贡献明细")
            lines.append("")
            lines.append("| 基因 | 位点 | 类别 | 贡献分 | 依据 |")
            lines.append("|------|------|------|--------|------|")
            for d in sorted(rare_details, key=lambda x: x.get("contribution", 0), reverse=True):
                notes = ", ".join(d.get("notes", []))
                var = d.get("variant", "")
                lines.append(
                    f"| {d.get('gene', '')} | `{var}` | 罕见功能变异 | {d.get('contribution', 0)} | {notes} |"
                )
            for h in gwas_hits:
                lines.append(
                    f"| {h.get('gene', '')} | {h.get('chrom', '')}:{h.get('pos', '')} "
                    f"({h.get('rsid', '')}) | GWAS lead SNP | +{h.get('points', 0)} | 基因型 {h.get('gt', '')} |"
                )
            lines.append("")
            lines.append(
                "> 贡献分为模型估算值，反映该变异/位点对当前表型的可能遗传贡献，非临床验证结果。"
            )
            lines.append("")
    else:
        lines.append("## 6. 综合风险评分")
        total = score_result.get("total_score", 0)
        level = score_result.get("risk_level", "无明确风险")
        meaning = score_result.get("risk_meaning", "")
        lines.append(f"- **总分**：{total}/100")
        lines.append(f"- **风险等级**：**{level}**")
        lines.append(f"- **含义**：{meaning}")
        lines.append("")
        lines.append("**分数组成**：")
        lines.append(f"- Tier 1 贡献：{comp.get('tier1_score', 0)} / {40}")
        lines.append(f"- Tier 2 贡献：{comp.get('tier2_score', 0)} / {25}")
        lines.append(f"- 文献支持：{comp.get('literature_score', 0)} / {20}")
        lines.append(f"- 人群罕见度：{round(comp.get('rarity_score', 0), 1)} / {10}")
        lines.append(f"- 性别/年龄校正：{score_result.get('sex_age_bonus', 0)} / {5}")
        if comp.get("min_gnomad_af") is not None:
            lines.append(f"- 最低 gnomAD AF：{comp['min_gnomad_af']}")
        lines.append("")

    lines.append("## 7. 数据局限性声明")
    lines.append("- 本分析基于 germline SNV/Indel，未覆盖 CNV/SV/表观遗传变异。")
    if is_complex:
        lines.append(
            "- 遗传贡献度评分用于量化遗传因素对复杂表型的可能解释比例，不等同于患病风险或临床诊断。"
        )
        lines.append(
            "- 复杂表型（如高尿酸、痛风、代谢综合征）受生活方式、环境和多基因累加影响，遗传贡献度仅反映当前可检出的遗传证据。"
        )
    else:
        lines.append("- 风险评分综合 ClinVar、OMIM、HPO 及文献证据，仅供科研参考，不构成临床诊断。")
    if vcf_qc and vcf_qc.get("common_variants_filtered"):
        lines.append(
            "- **输入 VCF 疑似已过滤常见变异**：锚定位点检出率较低，GWAS 常见风险等位基因在 VCF 中未保留。"
            "在已基因分型 callset 中可推断为 ref/ref；如需更准确的复杂表型遗传贡献评估，建议使用未过滤的完整 VCF 重新运行。"
        )
    lines.append("- 文献检索范围限定高影响因子英文期刊，可能存在发表偏倚。")
    lines.append("- 最终临床解读需结合完整家族史、表型及实验验证。")
    lines.append("")

    text = "\n".join(lines)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(text)
    return output_path
