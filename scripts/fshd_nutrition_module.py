#!/usr/bin/env python3
"""FSHD nutrition genetics standalone module.

Reads a short-read germline VCF and the curated FSHD nutrition gene set JSON,
queries known nutrigenetic SNPs and ClinVar P/LP variants in nutrient-pathway
genes, and produces a per-nutrient Markdown interpretation report.

The module is intentionally simpler than the full GPA pipeline: it does not
run VEP or compute PRS scores. Instead it reports genotype-nutrient associations
with explicit evidence tiers and clinical-check recommendations.

Usage:
    python3 fshd_nutrition_module.py \
        --vcf /path/to/sample.vcf.gz \
        --output-dir /path/to/output \
        --gene-set-json /path/to/fshd_nutrition_gene_set.json \
        --methylation-json /path/to/fshd_methylation_weights.json \
        --ancestry EAS
"""

import argparse
import json
import os
import shlex
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

BCFTOOLS = "/opt/homebrew/bin/bcftools"
CLINVAR_BED = "/Users/zhaorongli/.workbuddy/skills/dgra-prefilter/refs/clinvar_pathogenic_GRCh38.bed"


def run(cmd, check=True, capture=True):
    """Run a shell command and return stdout."""
    if isinstance(cmd, str):
        cmd = shlex.split(cmd)
    result = subprocess.run(
        cmd,
        shell=False,
        capture_output=capture,
        text=True,
        check=check,
    )
    if capture:
        return result.stdout
    return None


def quote(path):
    """Return a shell-safe quoted path."""
    return shlex.quote(path)


def load_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def build_gene_bed(gene_set, out_dir):
    """Build a BED file of all nutrient gene coordinates and curated SNP positions."""
    bed_path = os.path.join(out_dir, "fshd_nutrition_genes.bed")
    rows = []
    seen_genes = set()
    seen_snps = set()
    for nutrient, info in gene_set["nutrients"].items():
        for gene_entry in info.get("genes", []):
            coords = gene_entry.get("coords")
            if not coords:
                continue
            chrom = coords.get("chrom")
            start = coords.get("start")
            end = coords.get("end")
            gene = gene_entry["gene"]
            if chrom and start is not None and end and gene not in seen_genes:
                rows.append(f"{chrom}\t{start}\t{end}\t{gene}")
                seen_genes.add(gene)
        # Curated SNPs often lie outside the gene body (e.g., upstream/regulatory).
        # Include them with a +/- 1 bp window so they survive region extraction.
        for snp in info.get("gwas_lead_snps", []):
            chrom = snp.get("chrom")
            pos = snp.get("pos")
            if chrom and pos:
                start = pos - 1
                end = pos
                key = (chrom, start, end)
                if key not in seen_snps:
                    rows.append(f"{chrom}\t{start}\t{end}\t{snp.get('gene', '')}_SNP")
                    seen_snps.add(key)
    with open(bed_path, "w", encoding="utf-8") as f:
        f.write("\n".join(rows) + "\n")
    return bed_path


def extract_vcf_regions(vcf_path, bed_path, out_dir):
    """Extract sample variants falling within the nutrient gene regions."""
    out_vcf = os.path.join(out_dir, "fshd_nutrition_regions.vcf.gz")
    # Index output not needed for small region queries; just bgzip
    cmd = f"{BCFTOOLS} view -R {quote(bed_path)} -Oz -o {quote(out_vcf)} {quote(vcf_path)}"
    run(cmd)
    # Index for later -T queries
    run(f"{BCFTOOLS} index -t {quote(out_vcf)}")
    return out_vcf


def filter_clinvar_to_genes(gene_bed, out_dir):
    """Load ClinVar P/LP BED and keep only entries overlapping nutrient gene regions."""
    clinvar_filtered = os.path.join(out_dir, "fshd_nutrition_clinvar.bed")
    # Read gene intervals
    gene_intervals = defaultdict(list)
    with open(gene_bed, "r", encoding="utf-8") as f:
        for line in f:
            chrom, start, end, gene = line.strip().split("\t")
            gene_intervals[chrom].append((int(start), int(end), gene))

    # ClinVar BED is 0-based half-open: chrom, start, end, ...
    kept = []
    with open(CLINVAR_BED, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            chrom = parts[0]
            cstart = int(parts[1])
            cend = int(parts[2])
            # Check overlap with any gene on this chromosome
            for gstart, gend, gene in gene_intervals.get(chrom, []):
                if cstart < gend and cend > gstart:
                    kept.append(f"{chrom}\t{cstart}\t{cend}\t{gene}")
                    break
    with open(clinvar_filtered, "w", encoding="utf-8") as f:
        f.write("\n".join(kept) + "\n")
    return clinvar_filtered


def query_variants_in_bed(vcf_path, bed_path):
    """Return list of variants in the VCF that overlap the provided BED.

    BED must be 0-based half-open.
    """
    variants = []
    if os.path.getsize(bed_path) == 0:
        return variants
    cmd = f"{BCFTOOLS} view -H -R {quote(bed_path)} {quote(vcf_path)}"
    try:
        out = run(cmd)
    except subprocess.CalledProcessError:
        return variants
    for line in out.strip().split("\n"):
        if not line:
            continue
        cols = line.split("\t")
        chrom, pos, ref, alt = cols[0], cols[1], cols[3], cols[4]
        variants.append({
            "chrom": chrom,
            "pos": int(pos),
            "ref": ref,
            "alt": alt,
        })
    return variants


def query_clinvar_in_sample(vcf_path, clinvar_bed):
    """Query ClinVar P/LP variants found in the sample's nutrient gene regions."""
    # Build a 1-based inclusive position list for bcftools -T
    pos_bed = clinvar_bed.replace(".bed", "_1based.bed")
    with open(clinvar_bed, "r", encoding="utf-8") as f, open(pos_bed, "w", encoding="utf-8") as out:
        for line in f:
            chrom, start, end, gene = line.strip().split("\t")
            # 0-based [start, end) -> 1-based [start+1, end]
            out.write(f"{chrom}\t{int(start)+1}\t{end}\t{gene}\n")

    hits = []
    cmd = f"{BCFTOOLS} view -H -T {quote(pos_bed)} {quote(vcf_path)}"
    try:
        out = run(cmd)
    except subprocess.CalledProcessError:
        return hits
    for line in out.strip().split("\n"):
        if not line:
            continue
        cols = line.split("\t")
        chrom, pos, ref, alt = cols[0], cols[1], cols[3], cols[4]
        format_col = cols[8]
        sample_col = cols[9]
        gt_idx = format_col.split(":").index("GT")
        gt = sample_col.split(":")[gt_idx]
        # Determine zygosity (homozygous non-reference)
        alleles = gt.replace("|", "/").split("/")
        is_homozygous = (
            len(alleles) == 2
            and all(a != "0" and a != "." for a in alleles)
            and alleles[0] == alleles[1]
        )
        # Determine which gene(s) this position belongs to
        genes = []
        with open(clinvar_bed, "r", encoding="utf-8") as f:
            for l in f:
                c, s, e, g = l.strip().split("\t")
                if c == chrom and int(pos) - 1 >= int(s) and int(pos) <= int(e):
                    genes.append(g)
        hits.append({
            "chrom": chrom,
            "pos": int(pos),
            "ref": ref,
            "alt": alt,
            "genes": list(set(genes)),
            "gt": gt,
            "is_homozygous": is_homozygous,
        })
    return hits


def query_known_snp(vcf_path, snp):
    """Query a single known SNP position and return genotype info."""
    chrom = snp.get("chrom")
    pos = snp.get("pos")
    if not chrom or not pos:
        return None
    region = f"{chrom}:{pos}-{pos}"
    cmd = f"{BCFTOOLS} view -H -r {region} {quote(vcf_path)}"
    try:
        out = run(cmd)
    except subprocess.CalledProcessError:
        return None
    lines = [l for l in out.strip().split("\n") if l]
    if not lines:
        # Variant absent in VCF -> treat as ref/ref (0/0)
        return {
            "present": False,
            "gt": "0/0",
            "dosage": 0,
            "alleles": [],
        }
    line = lines[0]
    cols = line.split("\t")
    ref, alt = cols[3], cols[4]
    format_col = cols[8]
    sample_col = cols[9]
    gt_idx = format_col.split(":").index("GT")
    gt = sample_col.split(":")[gt_idx]
    # Dosage of effect allele
    effect = snp.get("effect_allele", "")
    alt_alleles = alt.split(",")
    alleles = [ref] + alt_alleles
    dosage = 0
    for allele in gt.replace("|", "/").split("/"):
        if allele == ".":
            continue
        if allele == "0":
            a = ref
        else:
            a = alt_alleles[int(allele) - 1]
        if a == effect:
            dosage += 1
    return {
        "present": True,
        "gt": gt,
        "dosage": dosage,
        "ref": ref,
        "alt": alt,
        "alleles": alleles,
    }


def methylation_modifier(methylation_data):
    """Return a label and modifier based on D4Z4 methylation."""
    if not methylation_data:
        return None, 1.0
    sample = methylation_data.get("sample_context", {})
    mean = sample.get("reported_mean_methylation_percent", 13.05)
    rules = methylation_data.get("modifier_rules", {}).get("overall_mean_methylation", {})
    if mean < 10:
        return rules.get("lt_10", {}).get("label", "very_low"), 1.2
    elif mean < 25:
        return rules.get("10_to_25", {}).get("label", "low"), 1.0
    elif mean < 40:
        return rules.get("25_to_40", {}).get("label", "intermediate"), 0.8
    else:
        return rules.get("ge_40", {}).get("label", "high"), 0.6


def interpret_nutrient(nutrient_key, nutrient_info, snp_results, clinvar_hits, modifier, ancestry="EAS"):
    """Generate a per-nutrient interpretation.

    Conservative interpretation rules:
    - "需关注" only when a curated known_pathogenic SNP (e.g., HFE C282Y/H63D)
      is present with the effect allele.
    - ClinVar P/LP variants are reported as supplementary findings and do NOT
      trigger the main nutrient status by themselves, because the raw ClinVar
      BED lacks phenotype, frequency, and review-status annotations.
    - Common GWAS risk SNPs contribute to "风险轻度/可能增加" based on
      effect-allele dosage.
    """
    label = nutrient_info.get("label", nutrient_key)
    tier = nutrient_info.get("tier", "C")
    clinical_check = nutrient_info.get("clinical_check", "")
    genes = nutrient_info.get("genes", [])

    # ClinVar P/LP hits in this nutrient's genes (reported, not used for status)
    nutrient_clinvar = [h for h in clinvar_hits if any(g["gene"] in h["genes"] for g in genes)]

    # Known SNPs split into pathogenic vs risk
    pathogenic_dosage = 0
    risk_dosage = 0
    for rsid, result in snp_results.items():
        if not result or result["dosage"] == 0:
            continue
        # Find the original SNP metadata to check variant_class
        snp_meta = next(
            (s for s in nutrient_info.get("gwas_lead_snps", []) if s["rsid"] == rsid), None
        )
        if snp_meta and snp_meta.get("variant_class") == "known_pathogenic":
            pathogenic_dosage += result["dosage"]
        else:
            risk_dosage += result["dosage"]

    # Tier C nutrients are exploratory: never trigger a high-impact alert from
    # GWAS or ClinVar findings. Only a curated known_pathogenic SNP can elevate.
    if tier == "C" and not pathogenic_dosage:
        status = "未见明确风险（探索性）"
        detail = f"未检测到明确的{label}高风险或致病性变异；该营养素关联证据较弱，仅作参考。"
    elif pathogenic_dosage >= 1:
        status = "需关注（检测到已知致病性营养基因变异）"
        detail = f"在{label}相关基因中检测到已知致病性变异，建议临床评估。"
    elif risk_dosage >= 2:
        status = "风险可能增加"
        detail = f"检测到多个{label}风险等位基因，可能降低代谢效率或增加需求。"
    elif risk_dosage == 1:
        status = "风险轻度增加"
        detail = f"检测到1个{label}风险等位基因，建议关注。"
    else:
        status = "未见明确风险"
        detail = f"未检测到已知的{label}高风险变异。"

    # Apply methylation modifier only for core muscle-skeletal nutrients
    if modifier != 1.0 and nutrient_key in {"vitamin_d", "calcium", "creatine", "protein_metabolism"}:
        if modifier > 1.0:
            detail += " 结合D4Z4低甲基化/严重肌肉受累背景，FSHD相关肌肉骨骼支持需求可能进一步增加。"
        elif modifier < 1.0:
            detail += " D4Z4甲基化相对不低，可作为肌肉受累严重程度的参考。"

    return {
        "label": label,
        "tier": tier,
        "status": status,
        "detail": detail,
        "clinical_check": clinical_check,
        "clinvar_hits": nutrient_clinvar,
        "snp_results": snp_results,
    }


def generate_markdown(report, sample_id, ancestry, methylation_label):
    """Generate the final Markdown report."""
    lines = [
        "# FSHD 肌肉骨骼营养遗传学报告",
        "",
        f"**样本**: {sample_id}",
        f"**祖先背景**: {ancestry}",
        f"**D4Z4甲基化背景**: {methylation_label or '未提供'}",
        "",
        "> **声明**：本报告为基因型-营养关联分析，不是临床诊断或营养处方。所有建议需由医生/注册营养师结合血清指标和临床评估确认。",
        "",
        "## 总体解读",
        "",
        "本报告基于FSHD患者肌肉-骨骼营养需求相关基因 curated list，分析维生素D、钙、维生素B12、叶酸、铁、肌酸等核心营养素，以及omega-3、蛋白质代谢、辅酶Q10等探索性营养素。",
        "",
        "| 营养素 | 证据等级 | 当前解读 | 建议检查/随访 |",
        "|--------|----------|----------|---------------|",
    ]
    for nutrient_key, interp in report.items():
        lines.append(
            f"| {interp['label']} | {interp['tier']} | {interp['status']} | {interp['clinical_check']} |"
        )

    lines.append("")
    lines.append("## 详细结果")
    lines.append("")

    for nutrient_key, interp in report.items():
        lines.append(f"### {interp['label']}（证据等级：{interp['tier']}）")
        lines.append("")
        lines.append(f"**解读**: {interp['detail']}")
        lines.append("")
        if interp['clinical_check']:
            lines.append(f"**建议检查**: {interp['clinical_check']}")
            lines.append("")
        if interp['snp_results']:
            lines.append("**已知SNP结果**:")
            lines.append("")
            lines.append("| rsID | 基因 | 效应等位基因 | 剂量 | 基因型 |")
            lines.append("|------|------|--------------|------|--------|")
            for rsid, res in interp['snp_results'].items():
                if res is None:
                    lines.append(f"| {rsid} | - | - | - | 坐标未解析 |")
                else:
                    lines.append(
                        f"| {rsid} | {res.get('gene', '-')} | "
                        f"{res.get('effect', '-')} | "
                        f"{res['dosage']} | {res['gt']} |"
                    )
            lines.append("")
        if interp['clinvar_hits']:
            lines.append("**ClinVar P/LP 附加发现**:")
            lines.append("")
            lines.append("| 位置 | 基因 | 基因型 | 纯合 | 备注 |")
            lines.append("|------|------|--------|------|------|")
            for hit in interp['clinvar_hits']:
                zyg = "纯合" if hit.get('is_homozygous') else "杂合"
                note = "ClinVar P/LP附加发现；与当前营养表型的相关性需临床确认"
                lines.append(
                    f"| {hit['chrom']}:{hit['pos']} {hit['ref']}>{hit['alt']} | "
                    f"{', '.join(hit['genes'])} | {hit.get('gt', '-')} | {zyg} | {note} |"
                )
            lines.append("")
        lines.append("")

    lines.append("## 限制与说明")
    lines.append("")
    lines.append("1. 常见变异效应量小，不能单独决定营养需求；必须结合血清检测和饮食评估。")
    lines.append("2. 本报告未纳入目前未知或尚未验证的营养-基因关联。")
    lines.append("3. 探索性营养素（Tier C，omega-3/蛋白质代谢/CoQ10）仅作参考，不用于临床决策。")
    lines.append("4. 东亚人群参考频率用于解释；若个体存在混合祖先，部分关联置信度降低。")
    lines.append("5. ClinVar P/LP附加发现未经表型、频率和审查星级过滤，仅作提示，不能单独作为营养诊断依据。")
    lines.append("")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="FSHD nutrition genetics analysis module")
    parser.add_argument("--vcf", required=True, help="Path to input VCF.gz")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--gene-set-json", required=True, help="Path to fshd_nutrition_gene_set.json")
    parser.add_argument("--methylation-json", default=None, help="Path to fshd_methylation_weights.json")
    parser.add_argument("--ancestry", default="EAS", help="Ancestry code for reporting")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    gene_set = load_json(args.gene_set_json)
    methylation_data = load_json(args.methylation_json) if args.methylation_json else None
    methylation_label, modifier = methylation_modifier(methylation_data)

    # Determine sample ID from VCF header
    sample_id = "unknown"
    header = run(f"{BCFTOOLS} view -h {quote(args.vcf)}")
    for line in header.strip().split("\n"):
        if line.startswith("#CHROM"):
            sample_id = line.split("\t")[-1]
            break

    # Build BED and extract regions
    gene_bed = build_gene_bed(gene_set, args.output_dir)
    filtered_vcf = extract_vcf_regions(args.vcf, gene_bed, args.output_dir)

    # Filter ClinVar to gene regions and query sample
    clinvar_bed = filter_clinvar_to_genes(gene_bed, args.output_dir)
    clinvar_hits = query_clinvar_in_sample(filtered_vcf, clinvar_bed)

    # Per-nutrient analysis
    report = {}
    for nutrient_key, nutrient_info in gene_set["nutrients"].items():
        snp_results = {}
        for snp in nutrient_info.get("gwas_lead_snps", []):
            rsid = snp["rsid"]
            res = query_known_snp(filtered_vcf, snp)
            if res is not None:
                res["gene"] = snp.get("gene", "")
                res["effect"] = snp.get("effect_allele", "")
            snp_results[rsid] = res

        report[nutrient_key] = interpret_nutrient(
            nutrient_key, nutrient_info, snp_results, clinvar_hits, modifier, args.ancestry
        )

    # Write report
    md = generate_markdown(report, sample_id, args.ancestry, methylation_label)
    report_path = os.path.join(args.output_dir, "fshd_nutrition_report.md")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(md)

    print(f"Report written to: {report_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
