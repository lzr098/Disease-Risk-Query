"""Markdown report generation for disease risk query.

This module generates a structured, disease-contribution-focused report.
It replaces the previous tier-centric layout with layered findings:
high-penetrance / moderate-penetrance / dosage-risk / GWAS-PRS / regulatory.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Any, Optional


LEVEL_MEANING: dict[str, str] = {
    "high": "高遗传贡献：当前可检出的遗传证据对该疾病解释力高",
    "moderate": "中等遗传贡献：当前可检出的遗传证据对该疾病有一定解释力",
    "low": "低遗传贡献：当前可检出的遗传证据有限，但存在潜在信号",
    "very_low": "当前可检出的遗传证据对该疾病解释力很低",
    "uncertain": "证据不足或无法评估",
    "none": "该层级未发现相关变异",
}


_LAYER_LABELS = {
    "mendelian_high": "高外显致病突变",
    "mendelian_mod": "中等外显变异",
    "known_pathogenic": "已知致病位点（模板）",
    "dosage_risk": "剂量风险位点",
    "gwas_prs": "GWAS / PRS 贡献",
    "regulatory": "调控区罕见变异",
    "clinvar_enriched": "ClinVar 注释变异",
}


def _get(v: dict, *keys: str) -> Any:
    """Get first matching key from variant dict (handles GPA lowercase output)."""
    for k in keys:
        if k in v:
            return v[k]
    return ""


def _variant_id(v: dict) -> str:
    chrom = _get(v, "CHROM", "chrom")
    pos = _get(v, "POS", "pos")
    ref = _get(v, "REF", "ref")
    alt = _get(v, "ALT", "alt")
    if chrom and pos and ref and alt:
        return f"{chrom}:{pos}:{ref}>{alt}"
    return v.get("variant") or ""


def _format_gt_status(inferred: bool, gt: str = "") -> str:
    if inferred:
        return "推断 ref/ref"
    if gt in ("0/0", "0|0"):
        return "样本检出 ref/ref"
    return "样本检出"


def _variant_table_row(v: dict, show_clinvar: bool = False, gene_context: dict[str, str] | None = None) -> str:
    # v may be a contribution wrapper with "raw" pointing to the GPA variant,
    # or a raw GPA variant dict.
    raw = v.get("raw", v)
    gene = _get(raw, "GENE", "gene")
    var_id = _variant_id(raw) or v.get("variant", "")
    hgvsp = _get(raw, "HGVSp", "hgvsp", "primary_hgvsp")
    hgvsc = _get(raw, "HGVSc", "hgvsc", "primary_hgvsc")
    clin_sig = _get(raw, "clinvar_sig", "CLIN_SIG", "clinvar")
    impact = _get(raw, "IMPACT", "impact", "primary_impact")
    cons = _get(raw, "Consequence", "consequence", "primary_consequence")
    gt = _get(raw, "GT", "gt", "zygosity")
    gnomad = _get(raw, "gnomAD_AF", "gnomad_af")
    flags = ", ".join(raw.get("_drq_flags", [])) or "-"

    gene_contrib = raw.get("_gene_contribution")
    gene_pen = raw.get("_gene_penetrance")
    contrib_str = f"{gene_contrib}" if gene_contrib is not None else "-"
    pen_str = gene_pen or "-"

    ctx = ""
    if gene_context and gene in gene_context:
        ctx = gene_context[gene]

    if show_clinvar:
        return (
            f"| {gene} | {ctx} | `{var_id}` | {hgvsc} | {hgvsp} | "
            f"{clin_sig} | {impact}/{cons} | {gt} | {gnomad} | {contrib_str} | {pen_str} | {flags} |"
        )
    return (
        f"| {gene} | {ctx} | `{var_id}` | {hgvsc} | {hgvsp} | "
        f"{impact}/{cons} | {gt} | {gnomad} | {contrib_str} | {pen_str} | {flags} |"
    )


def _known_variant_table_row(d: dict) -> str:
    rsid = d.get("rsid") or "-"
    gene = d.get("gene") or "-"
    variant = d.get("variant") or "-"
    risk_allele = d.get("risk_allele") or d.get("effect_allele") or "-"
    dosage = d.get("dosage") if d.get("dosage") is not None else "-"
    gt = d.get("gt") or "0/0"
    status = _format_gt_status(d.get("inferred_ref_ref", False), gt)
    or_val = d.get("or_per_allele") or d.get("or_value")
    or_str = f"{or_val:.2f}" if isinstance(or_val, (int, float)) else "-"
    beta = d.get("beta")
    beta_str = f"{beta:.4f}" if isinstance(beta, (int, float)) else "-"
    contrib = d.get("contribution")
    contrib_str = f"{contrib:.3f}" if isinstance(contrib, (int, float)) else "-"
    note = d.get("note", "")
    return (
        f"| {rsid} | {gene} | `{variant}` | {gt} | {status} | {risk_allele} | "
        f"{dosage} | {or_str} | {beta_str} | {contrib_str} | {note} |"
    )


def _render_layer_table(layer: str, items: list[dict], gene_context: dict[str, str] | None = None) -> list[str]:
    lines: list[str] = []
    if not items:
        lines.append("_未检出_")
        return lines
    ctx = gene_context or {}

    if layer == "gwas_prs":
        lines.append("| SNP | 基因 | 位点 | 基因型 | 来源 | 风险等位 | 剂量 | OR | beta | 贡献 | 说明 |")
        lines.append("|-----|------|------|--------|------|----------|------|----|----|------|------|")
        for d in items:
            lines.append(_known_variant_table_row(d))
    elif layer in ("dosage_risk", "known_pathogenic"):
        lines.append("| SNP | 基因 | 功能 | 位点 | 基因型 | 来源 | 风险等位 | 剂量 | 贡献 | 置信度 | 说明 |")
        lines.append("|-----|------|------|------|--------|------|----------|------|------|--------|------|")
        for d in items:
            rsid = d.get("rsid") or "-"
            gene = d.get("gene") or "-"
            # Get gene context from known variant note or gene_context
            gene_fn = ctx.get(gene, "")
            if not gene_fn and d.get("note"):
                # Extract first sentence of note for function context
                note = d.get("note", "")
                gene_fn = note.split("；")[0] if "；" in note else note[:60]
            variant = d.get("variant") or "-"
            risk_allele = d.get("risk_allele") or d.get("effect_allele") or "-"
            dosage = d.get("dosage") if d.get("dosage") is not None else "-"
            gt = d.get("gt") or "0/0"
            status = "推断 ref/ref" if d.get("inferred_ref_ref") else "样本检出"
            contrib = d.get("contribution")
            contrib_str = f"{contrib:.3f}" if isinstance(contrib, (int, float)) else "-"
            confidence = d.get("confidence", "-")
            note_line = d.get("note", "")
            lines.append(
                f"| {rsid} | {gene} | {gene_fn} | `{variant}` | {gt} | {status} | {risk_allele} | "
                f"{dosage} | {contrib_str} | {confidence} | {note_line} |"
            )
    else:
        lines.append("| 基因 | 基因-表型关系 | 位点 | cDNA | 蛋白 | ClinVar | 影响/类型 | 合子性 | gnomAD AF | 基因贡献 | 外显率 | 标志 |")
        lines.append("|------|--------------|------|------|------|---------|-----------|--------|-----------|----------|--------|------|")
        show_clinvar = layer in ("mendelian_high", "mendelian_mod", "clinvar_enriched")
        for v in items:
            lines.append(_variant_table_row(v, show_clinvar=show_clinvar, gene_context=ctx))
    return lines


def _executive_summary(
    score_result: dict,
    contribution: dict,
    known_genotypes: list[dict],
) -> list[str]:
    lines: list[str] = []

    # Key findings
    findings: list[str] = []
    high = contribution.get("mendelian_high", [])
    mod = contribution.get("mendelian_mod", [])
    known = contribution.get("known_pathogenic", [])
    dosage = contribution.get("dosage_risk", [])
    gwas = contribution.get("gwas_prs", {})

    if high:
        genes = sorted({h.get("gene", "") for h in high})
        findings.append(f"发现高外显致病突变（{', '.join(genes)}）")
    elif known:
        real_known = [k for k in known if not k.get("inferred_ref_ref")]
        if real_known:
            detail_parts: list[str] = []
            for k in real_known[:3]:
                rsid = k.get("rsid", "")
                gene = k.get("gene", "")
                gt = k.get("gt", "")
                note = k.get("note", "")
                # Build concise description
                desc = f"{gene} {rsid}"
                if gt:
                    desc += f"（{gt}）"
                if note:
                    # Trim note to keep it concise
                    note_short = note[:80] + ("..." if len(note) > 80 else "")
                    desc += f"：{note_short}"
                detail_parts.append(desc)
            detail = "；".join(detail_parts)
            findings.append(f"发现已知致病位点真实检出：{detail}")
    elif mod:
        genes = sorted({m.get("gene", "") for m in mod})
        findings.append(f"发现中等外显变异（{', '.join(genes)}）")

    real_dosage = [d for d in dosage if not d.get("inferred_ref_ref")]
    if real_dosage:
        items = [f"{d.get('rsid')} (dosage={d.get('dosage')})" for d in real_dosage[:3]]
        findings.append(f"剂量风险位点真实检出：{', '.join(items)}")

    real_gwas = [v for v in gwas.get("variants", []) if not v.get("inferred_ref_ref")]
    if real_gwas:
        findings.append(f"GWAS/PRS 位点真实检出 {len(real_gwas)} 个")

    if not findings:
        findings.append("本次分析未发现明确的高外显致病突变或强风险等位基因")

    # Mention inferred ref/ref count
    inferred_known = [kg for kg in known_genotypes if kg.get("sample_gt", {}).get("inferred_ref_ref")]
    if inferred_known:
        findings.append(
            f"{len(inferred_known)} 个已知风险位点未在 VCF 中保留，已按 ref/ref（0/0）推断"
        )

    for f in findings[:5]:
        lines.append(f"- {f}")

    lines.append("")
    lines.append(
        "> **注意**：未检出高外显致病突变不等于排除该疾病的遗传病因。"
        "本报告仅反映当前 VCF 可检出的遗传证据。"
    )
    return lines


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
    apoe_result: Optional[dict] = None,  # deprecated, kept for compatibility
    gwas_summary: Optional[dict] = None,
    literature_summary: Optional[dict] = None,
    disease_reference: Optional[dict] = None,
    gwas_lead_snps: Optional[list[dict]] = None,
    vcf_qc: Optional[dict] = None,
    disease_mode: str = "mendelian",
    domain_dive_candidates: Optional[list[dict]] = None,
    disease_space: Optional[dict] = None,
) -> Path:
    """Generate the final Markdown disease contribution report.

    Args:
        disease_space: Optional dict with variant analysis counts:
            {'total_variants': 5366066, 'analyzed_variants': 6229,
             'known_variants_queried': 48, 'known_found': 17,
             'known_inferred': 31}
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    contribution = score_result.get("contribution", score_result)
    layer_levels = contribution.get("layer_levels", {})

    lines: list[str] = []
    total = score_result.get("total_score", 0)
    level = score_result.get("overall_level", "uncertain")
    meaning = LEVEL_MEANING.get(level, level)

    # --- Header with score and breakdown ---
    lines.append(f"# 疾病遗传贡献度评估报告：{disease_name} — **{total}/100** ({level})")
    lines.append("")
    lines.append(f"**{total}/100 — {level}**：{meaning}")
    lines.append("")

    # Score composition (weights from contribution_scorer.py, not in result dict)
    layer_scores = {}
    for layer in ("mendelian_high", "mendelian_mod", "known_pathogenic", "dosage_risk", "gwas_prs", "regulatory"):
        score = sum(
            x.get("contribution", 0) for x in contribution.get(layer, [])
            if isinstance(x, dict)
        )
        if layer == "gwas_prs":
            score = abs(contribution.get("gwas_prs", {}).get("score", 0.0))
        layer_scores[layer] = score

    parts: list[str] = []
    for layer in ("mendelian_high", "mendelian_mod", "known_pathogenic", "dosage_risk", "gwas_prs", "regulatory"):
        s = layer_scores[layer]
        if s > 0 or layer == "known_pathogenic":
            label = _LAYER_LABELS.get(layer, layer)
            parts.append(f"{label}: {s:.3f}")
    if parts:
        lines.append(f"**分数组成**：{' | '.join(parts)}")
        lines.append(f"**总分** = Σ min(layer, 1.0) × weight = {contribution.get('overall_score', total/100):.3f}")
    lines.append("")

    lines.append("## 1. 执行摘要")
    lines.extend(_executive_summary(score_result, contribution, gwas_lead_snps or []))
    lines.append("")

    # 2. Query summary
    lines.append("## 2. 查询与样本信息")
    lines.append(f"- **目标疾病**：{disease_name}")
    lines.append(f"- **评估模式**：{disease_mode}")
    lines.append(f"- **HPO ID**：{hpo_id or 'N/A'} ({hpo_name or 'N/A'})")
    lines.append(f"- **样本性别**：{sex} | **年龄**：{age if age is not None else '未知'}")
    lines.append(f"- **报告生成时间**：{datetime.now().isoformat(timespec='minutes')}")
    lines.append("")

    lines.append("### VCF 质量与分析统计")
    ds = disease_space or {}
    if ds:
        total_vcf = ds.get("total_variants", vcf_qc.get("total_variants", "N/A") if vcf_qc else "N/A")
        analyzed = ds.get("analyzed_variants", "N/A")
        kv_queried = ds.get("known_variants_queried", 0)
        kv_found = ds.get("known_found", 0)
        kv_inferred = ds.get("known_inferred", 0)
        lines.append(f"- **VCF 总变异数**：{total_vcf:,} 个")
        lines.append(f"- **本次分析纳入**：{analyzed:,} 个（过滤至疾病相关基因区域）")
        if kv_queried:
            lines.append(f"- **关键 SNP 查询**：{kv_queried} 个已知风险位点，其中 {kv_found} 个真实检出 + {kv_inferred} 个未 call 到（按 ref/ref 推断）")
        lines.append(
            "- **解读原则**：未 call 到的位点一律推断为 ref/ref（0/0），"
            "需注意前序 genotyping 及过滤流程的可靠性。"
            "常见 SNP 未保留提示上游可能经过硬过滤，多基因贡献评估可能存在低估。"
        )
    elif vcf_qc and vcf_qc.get("checked"):
        lines.append(
            f"- **锚定位点真实检出率**：{vcf_qc.get('anchor_snps_present', 0)} / "
            f"{vcf_qc.get('anchor_snps_checked', 0)} "
            f"({vcf_qc.get('presence_rate', 0):.0%})"
        )
        lines.append(f"- **VCF 总变异数（近似）**：{vcf_qc.get('total_variants', 'N/A'):,}")
        if vcf_qc.get("common_variants_filtered") or vcf_qc.get("likely_filtered"):
            lines.append(
                "- **常见变异过滤提示**：该 VCF 的锚定位点真实检出率低于阈值，"
                "说明常见 SNP 位点未在 VCF 中保留。对于已基因分型的 callset，"
                "未保留的位点一律推断为 ref/ref（0/0）；GWAS/PRS 维度仍按完整权重评估，但可信度下降。"
            )
        else:
            lines.append("- **VCF 完整性**：锚定位点真实检出率正常，GWAS/PRS 维度可正常评估。")
    else:
        lines.append(f"- **VCF QC**：{vcf_qc.get('note', '未执行') if vcf_qc else '未执行'}")
    lines.append("")

    # 3. Disease profile summary
    lines.append("## 3. 疾病模板摘要")
    if disease_reference:
        meta = disease_reference.get("metadata", {})
        counts = meta.get("counts", {})
        lines.append(f"- **核心基因**：{counts.get('core_genes', 0)} 个")
        lines.append(f"- **已知 ClinVar 致病突变**：{counts.get('clinvar_variants', 0)} 个")
        lines.append(f"- **GWAS lead SNPs**：{counts.get('gwas_snps', 0)} 个")
        lines.append(f"- **关键文献条目**：{counts.get('literature_entries', 0)} 条")
    sources = gene_set_result.get("sources", {})
    lines.append(f"- **基因来源**：OMIM {sources.get('omim', 0)} 个，HPO {sources.get('hpo', 0)} 个，文献扩展 {sources.get('extra', 0)} 个")
    lines.append(f"- **合并去重后基因总数**：{gene_set_result.get('total', 0)}")
    genes = gene_set_result.get("merged_genes", [])
    if genes:
        lines.append(f"- **基因列表（前 30）**：{', '.join(genes[:30])}{' ...' if len(genes) > 30 else ''}")
    lines.append("")

    # Enrich variants with template-based gene contribution / penetrance
    gene_contribution_map = (disease_reference or {}).get("gene_contribution_map", {})
    gene_penetrance_map = (disease_reference or {}).get("gene_penetrance_map", {})
    for v in (
        gpa_result.get("tier1_variants", [])
        + gpa_result.get("tier2_variants", [])
        + gpa_result.get("tier3_variants", [])
    ):
        gene = _get(v, "GENE", "gene")
        if gene and gene in gene_contribution_map:
            v["_gene_contribution"] = gene_contribution_map[gene]
            v["_gene_penetrance"] = gene_penetrance_map.get(gene, "")

    # 4. Layered findings
    lines.append("## 4. 分层发现")
    lines.append("")

    layer_order = [
        "mendelian_high", "mendelian_mod", "known_pathogenic",
        "dosage_risk", "gwas_prs", "regulatory",
    ]
    raw_layer_items = {
        "mendelian_high": contribution.get("mendelian_high", []),
        "mendelian_mod": contribution.get("mendelian_mod", []),
        "known_pathogenic": contribution.get("known_pathogenic", []),
        "dosage_risk": contribution.get("dosage_risk", []),
        "gwas_prs": contribution.get("gwas_prs", {}).get("variants", []),
        "regulatory": contribution.get("regulatory", []),
    }
    # For GWAS/PRS, suppress all-ref/ref zero-contribution noise in the table,
    # but keep the summary counts honest.
    # For known_pathogenic, also filter out absent (ref/ref) to reduce noise.
    layer_items: dict[str, list[dict]] = {}
    for layer in layer_order:
        items = list(raw_layer_items[layer])
        if layer == "gwas_prs":
            items = [
                x for x in items
                if not x.get("inferred_ref_ref") or x.get("contribution", 0) != 0
            ]
        if layer == "known_pathogenic":
            items = [
                x for x in items
                if not x.get("inferred_ref_ref") or x.get("contribution", 0) != 0
            ]
        layer_items[layer] = items

    layer_scores = {
        "mendelian_high": sum(x.get("contribution", 0) for x in raw_layer_items["mendelian_high"]),
        "mendelian_mod": sum(x.get("contribution", 0) for x in raw_layer_items["mendelian_mod"]),
        "known_pathogenic": sum(x.get("contribution", 0) for x in raw_layer_items["known_pathogenic"]),
        "dosage_risk": sum(x.get("contribution", 0) for x in raw_layer_items["dosage_risk"]),
        "gwas_prs": abs(contribution.get("gwas_prs", {}).get("score", 0.0)),
        "regulatory": sum(x.get("contribution", 0) for x in raw_layer_items["regulatory"]),
    }

    section_counter = 1

    # Build gene-context lookup from disease reference and known variant notes
    gene_context: dict[str, str] = {}
    if disease_reference:
        core_genes_set = set(disease_reference.get("core_genes", []))
        for g in core_genes_set:
            gene_context[g] = "核心基因"
    # Enrich with known variant notes (first sentence as function context)
    for k in contribution.get("known_pathogenic", []):
        gene = k.get("gene", "")
        note = k.get("note", "")
        if gene and note and gene not in gene_context:
            gene_context[gene] = note.split("；")[0] if "；" in note else note[:60]

    for layer in layer_order:
        label = _LAYER_LABELS[layer]
        items = layer_items[layer]
        raw_items = raw_layer_items[layer]
        score = layer_scores[layer]
        level = layer_levels.get(layer, "none")
        real_count = sum(1 for x in raw_items if not x.get("inferred_ref_ref"))
        inferred_count = len(raw_items) - real_count

        lines.append(f"### 4.{section_counter} {label}")
        section_counter += 1
        lines.append(f"- **层级得分**：{score:.3f}")
        lines.append(f"- **层级等级**：{level}（{LEVEL_MEANING.get(level, level)}）")
        lines.append(f"- **真实检出**：{real_count} 个 | **推断 ref/ref**：{inferred_count} 个")
        lines.append("")
        if layer in ("gwas_prs", "dosage_risk") and not items and raw_items:
            # All known risk variants are ref/ref with zero contribution; show compact summary
            label_short = "GWAS/PRS" if layer == "gwas_prs" else "剂量风险"
            lines.append(
                f"_该疾病共定义 {len(raw_items)} 个 {label_short} 位点，"
                f"本次样本中 {inferred_count} 个未真实检出（已推断为 ref/ref），"
                "无风险等位基因贡献。_"
            )
        else:
            lines.extend(_render_layer_table(layer, items, gene_context))
        lines.append("")

    # Other notable Tier 2/3 variants
    section_num = 5
    lines.append(f"## {section_num}. 其他值得关注的 Tier 2/3 变异")
    section_num += 1
    tier1 = gpa_result.get("tier1_variants", [])
    tier2 = gpa_result.get("tier2_variants", [])
    tier3 = gpa_result.get("tier3_variants", [])

    notable: list[dict] = []
    for v in tier2 + tier3:
        if v.get("clinvar_plp"):
            notable.append(v)
            continue
        if _get(v, "IMPACT", "impact", "primary_impact").upper() == "HIGH":
            notable.append(v)
            continue
        gene = _get(v, "GENE", "gene")
        if gene and gene in (gene_contribution_map or {}):
            notable.append(v)

    # Remove duplicates already shown in layered findings
    shown_keys = {_variant_id(x) for x in tier1}
    shown_keys.update({_variant_id(x) for x in layer_items["mendelian_mod"]})
    shown_keys.update({_variant_id(x) for x in layer_items["regulatory"]})
    notable = [v for v in notable if _variant_id(v) not in shown_keys]

    if notable:
        # Build gene-context lookup from disease reference
        core_genes_set = set(disease_reference.get("core_genes", [])) if disease_reference else set()
        lines.append(f"以下 {len(notable)} 个变异位于疾病核心基因内，建议结合基因-表型关系关注：")
        lines.append("")
        lines.append("| 基因 | 与疾病关系 | 位点 | 蛋白 | ClinVar | 影响 | 合子性 | gnomAD AF |")
        lines.append("|------|------------|------|------|---------|------|--------|-----------|")
        for v in notable:
            gene = _get(v, "GENE", "gene")
            gene_ctx = "核心基因" if gene in core_genes_set else "疾病相关"
            var_id = _variant_id(v)
            hgvsp = _get(v, "HGVSp", "hgvsp", "primary_hgvsp")
            clin_sig = _get(v, "clinvar_sig", "CLIN_SIG", "clinvar")
            impact = _get(v, "IMPACT", "impact", "primary_impact")
            cons = _get(v, "Consequence", "consequence", "primary_consequence")
            gt = _get(v, "GT", "gt", "zygosity")
            gnomad = _get(v, "gnomAD_AF", "gnomad_af")
            lines.append(
                f"| {gene} | {gene_ctx} | `{var_id}` | {hgvsp} | "
                f"{clin_sig} | {impact}/{cons} | {gt} | {gnomad} |"
            )
        lines.append("")
    else:
        lines.append("_未检出其他需要特别关注的 Tier 2/3 变异_")
    lines.append("")

    # ClinVar enriched summary
    clinvar_enriched = contribution.get("clinvar_enriched", [])
    if clinvar_enriched:
        lines.append(f"## {section_num}. ClinVar 注释概览")
        section_num += 1
        lines.append(f"- 疾病模板基因区内共有 **{len(clinvar_enriched)}** 个变异带有 ClinVar 注释。")
        plp = [v for v in clinvar_enriched if v.get("clinvar_category") in ("pathogenic", "likely_pathogenic")]
        vus = [v for v in clinvar_enriched if v.get("clinvar_category") == "vus"]
        other = [v for v in clinvar_enriched if v.get("clinvar_category") not in ("pathogenic", "likely_pathogenic", "vus")]
        lines.append(f"- P/LP：{len(plp)} 个 | VUS：{len(vus)} 个 | 其他：{len(other)} 个")
        lines.append("- ClinVar 分类仅作为信息增强，VUS 等不自动降权。")
        lines.append("")

    # Domain-dive candidates
    if domain_dive_candidates:
        lines.append(f"## {section_num}. Domain-dive 升级候选分析")
        section_num += 1
        lines.append(f"- 对 Tier 2/3 中的核心基因变异进行了疾病特异性结构域深度分析，发现 **{len(domain_dive_candidates)} 个** 潜在升级候选。")
        lines.append("- 以下位点落在该疾病关键功能区域或邻近关键残基，建议人工复核。")
        lines.append("")
        lines.append("| 基因 | 蛋白变化 | 关键区域 | 最近关键残基 | 推荐 | 依据 |")
        lines.append("|------|----------|----------|--------------|------|------|")
        for item in domain_dive_candidates:
            v = item["variant"]
            dd = item["domain_dive"]
            gene = _get(v, "GENE", "gene")
            hgvsp = _get(v, "HGVSp", "hgvsp", "primary_hgvsp")
            regions = ", ".join(r["name"] for r in dd.get("matched_regions", [])) or "-"
            nearby = dd.get("nearest_critical_residues", [])
            nearby_str = "; ".join(f"{n['residue']}({n['distance']}aa)" for n in nearby[:3]) if nearby else "-"
            rec = dd.get("upgrade_recommendation", "")
            rec_zh = {"tier2_candidate": "考虑升 Tier 2", "monitor": "密切监测", "no_evidence": "无证据"}.get(rec, rec)
            reasoning = dd.get("reasoning", "")
            lines.append(f"| {gene} | {hgvsp} | {regions} | {nearby_str} | {rec_zh} | {reasoning} |")
        lines.append("")

    # Data limitations
    lines.append(f"## {section_num}. 数据局限性与建议")
    lines.append("- 本分析基于 germline SNV/Indel，未覆盖 CNV/SV/表观遗传变异。")
    lines.append(
        "- **未检出高外显致病突变不等于排除遗传病因**：部分致病变异可能位于当前 VCF "
        "未覆盖区域、未被注释的基因，或需要家系/功能实验验证。"
    )
    lines.append("- 遗传贡献度评分反映当前可检出的遗传证据，不等同于患病风险或临床诊断。")
    if vcf_qc and (vcf_qc.get("common_variants_filtered") or vcf_qc.get("likely_filtered")):
        lines.append(
            "- **输入 VCF 疑似已过滤常见变异**：锚定位点真实检出率较低，GWAS/PRS 常见风险等位基因部分未保留。"
            "推断 ref/ref 会在一定程度上低估多基因风险；如需更准确评估，建议使用未过滤的完整 VCF。"
        )
    lines.append("- 文献检索范围可能受数据库更新频率和发表偏倚影响。")
    lines.append("- 最终临床解读需结合完整家族史、表型及实验验证。")
    lines.append("")

    text = "\n".join(lines)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(text)
    return output_path
