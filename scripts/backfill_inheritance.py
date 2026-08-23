#!/usr/bin/env python3
"""给疾病模板的 mendelian 层基因自动补 inheritance（遗传模式）字段。

数据源：
  1. ClinGen gene-disease validity curations（逐记录，含 MOI + CLASSIFICATION）
     路径见 CLINGEN_CSV，本地文件，CSV 带引号格式，表头在含 'GENE SYMBOL' 的行。
  2. HPO genes_to_phenotype.txt（gene → inheritance 术语）

补字段规则（重要，记录于此）：
  R1. 范围：仅 mendelian_high / mendelian_mod 层基因。
      strong_gwas / gwas 层不补——复杂病/多基因病没有单基因遗传模式。
  R2. 数据源合并：ClinGen 的 MOI（经过可信度过滤）∪ HPO 的 inheritance 术语。
  R3. ClinGen 可信度过滤：排除 Disputed / Refuted / Contested /
      No Known Disease Relationship / Not Yet Evaluated 等否定或未定分类，
      仅保留 Definitive / Strong / Moderate / Limited。
  R4. 填值：
      - 合并后空集   → 不填（保持空，表示未知或复杂病背景）
      - 单值         → 填该值（如 "AR"）
      - 多值         → 按固定顺序逗号连接（如 "AD,AR"），不硬猜主模式。
        多值反映真实情况：同一基因通过不同疾病可对应多种遗传模式，
        下游评分器应对多值做保守处理（杂合时同时考虑 AD 致病 / AR 携带）。
  R5. 覆盖策略：统一用数据库值覆盖所有 mendelian 层基因的 inheritance
      （包括历史上手工填写的），保证标准一致；差异在 stdout 输出供审阅。
  R6. 只写回 gene_set 的 inheritance 字段，其余字段、注释头、TEMPLATE_NAME 原样保留。

用法：
  python scripts/backfill_inheritance.py [--dry-run]
"""
import csv
import importlib
import os
import pprint
import re
import sys
from collections import defaultdict

# 路径
SKILL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TPL_DIR = os.path.join(SKILL_DIR, "references", "disease_templates")
CLINGEN_CSV = os.path.expanduser(
    "~/.workbuddy/data/clingen/ClinGen_gene_validity_curations.tsv"
)
HPO_TXT = "/tmp/hpo_genes_to_phenotype.txt"

# ClinGen 可信度：保留 vs 排除
KEEP_CLASSES = {"Definitive", "Strong", "Moderate", "Limited"}
# HPO inheritance 术语 -> 标准缩写
HPO_INH = {
    "HP:0000006": "AD",   # Autosomal dominant
    "HP:0000007": "AR",   # Autosomal recessive
    "HP:0001417": "XL",   # X-linked
    "HP:0001423": "XLd",  # X-linked dominant
    "HP:0001419": "XLr",  # X-linked recessive
    "HP:0001427": "MT",   # Mitochondrial
    "HP:0001450": "YL",   # Y-linked
    "HP:0001428": "SD",   # Somatic
    "HP:0032113": "SD",   # Semidominant
}
# 多值时的固定输出顺序
MOI_ORDER = ["AD", "AR", "XL", "XLd", "XLr", "MT", "YL", "SD"]


def load_clingen():
    """逐记录读 ClinGen，返回 gene -> set(MOI)，仅保留可信分类。"""
    gene_moi = defaultdict(set)
    if not os.path.exists(CLINGEN_CSV):
        print(f"[WARN] ClinGen 文件不存在: {CLINGEN_CSV}", file=sys.stderr)
        return gene_moi
    with open(CLINGEN_CSV) as f:
        rd = csv.reader(f)
        for row in rd:
            if row and row[0] == "GENE SYMBOL":
                break
        for row in rd:
            if len(row) < 7:
                continue
            gene, moi, cls = row[0].strip(), row[4].strip(), row[6].strip()
            if not gene or not moi or moi in ("+++++++++", "UD", ""):
                continue
            if cls in KEEP_CLASSES:
                gene_moi[gene].add(moi)
    return gene_moi


def load_hpo():
    """读 HPO genes_to_phenotype.txt，返回 gene -> set(MOI)。"""
    gene_moi = defaultdict(set)
    if not os.path.exists(HPO_TXT):
        print(f"[WARN] HPO 文件不存在: {HPO_TXT}", file=sys.stderr)
        return gene_moi
    with open(HPO_TXT) as f:
        next(f)  # header
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            gene, hpo_id = p[1].strip(), p[2].strip()
            if hpo_id in HPO_INH:
                gene_moi[gene].add(HPO_INH[hpo_id])
    return gene_moi


def resolve(gene, clingen, hpo):
    """返回该基因的遗传模式（排序后的列表）。

    策略（记录）：ClinGen 优先。
      - 基因在 ClinGen 有可信记录 → 只用 ClinGen 的 MOI（curator 评审，权威）。
      - 无 ClinGen 记录 → 用 HPO 兜底。
    理由：HPO 的 genes_to_phenotype 存在噪声——同一 OMIM 疾病会被
    同时标注矛盾的遗传模式（如 AD 和 AR），若与 ClinGen 简单取并集会
    把噪声混入。ClinGen 有 curator 逐条评审 + 可信度分级，更可靠。
    """
    if gene in clingen:
        return sorted(clingen[gene], key=lambda m: MOI_ORDER.index(m) if m in MOI_ORDER else 99)
    return sorted(hpo.get(gene, set()), key=lambda m: MOI_ORDER.index(m) if m in MOI_ORDER else 99)


def fmt(mois):
    return ",".join(mois) if mois else ""


def rewrite_template(module_name, tempplate_obj):
    """保留注释头 + TEMPLATE_NAME 行，重新生成 TEMPLATE 并写回。"""
    path = os.path.join(TPL_DIR, module_name + ".py")
    with open(path) as f:
        content = f.read()
    # 精确定位 TEMPLATE_NAME 赋值行（避免误匹配 docstring 里的 "TEMPLATE_NAME" 字样）
    m = re.search(r'\nTEMPLATE_NAME\s*=\s*["\']', content)
    if not m:
        raise RuntimeError(f"{module_name}: 未找到 TEMPLATE_NAME 赋值行")
    # header = 赋值行之前的内容（docstring 等）+ 该赋值行本身（到行尾）
    line_end = content.index("\n", m.start() + 1) + 1
    header = content[:line_end]
    # 重建 TEMPLATE 部分（统一 pprint 格式）
    body = "TEMPLATE = " + pprint.pformat(tempplate_obj, width=120, sort_dicts=False) + "\n"
    with open(path, "w") as f:
        f.write(header + body)


def main():
    dry_run = "--dry-run" in sys.argv
    clingen = load_clingen()
    hpo = load_hpo()
    print(f"[数据源] ClinGen 可信 MOI 基因数: {len(clingen)}, HPO 基因数: {len(hpo)}")

    files = sorted(
        f[:-3] for f in os.listdir(TPL_DIR)
        if f.endswith(".py") and f != "__init__.py"
    )

    # 全局统计
    stat_single = stat_multi = stat_none = stat_skip = 0
    conflicts = []  # 手填 vs 自动 差异
    changed = []  # (模板, 基因, 旧值, 新值)

    for name in files:
        # 清缓存重新导入
        for m in list(sys.modules):
            if m == name:
                del sys.modules[m]
        sys.path.insert(0, TPL_DIR)
        mod = importlib.import_module(name)
        T = mod.TEMPLATE
        tpl_changed = []
        for g in T.get("gene_set", []):
            tier = g.get("tier", "")
            if tier not in ("mendelian_high", "mendelian_mod"):
                stat_skip += 1
                continue
            gene = g.get("gene", "")
            mois = resolve(gene, clingen, hpo)
            new_val = fmt(mois)
            old_val = g.get("inheritance", "")
            if not new_val:
                stat_none += 1
                # 无数据库记录：保持原状（可能已有手填）
                continue
            if len(mois) == 1:
                stat_single += 1
            else:
                stat_multi += 1
            if old_val != new_val:
                tpl_changed.append((gene, old_val, new_val))
                if old_val:  # 有手填但被覆盖，记录差异
                    conflicts.append((name, gene, old_val, new_val))
            g["inheritance"] = new_val

        if tpl_changed:
            changed.append((name, tpl_changed))
        if not dry_run:
            rewrite_template(name, T)
        print(f"  {name}: 更新 {len(tpl_changed)} 个基因")

    print()
    print(f"=== 统计 ===")
    print(f"  单值填充: {stat_single}")
    print(f"  多值填充(逗号连接): {stat_multi}")
    print(f"  无数据库记录(留空): {stat_none}")
    print(f"  跳过(非mendelian层): {stat_skip}")
    print(f"  总更新记录数: {sum(len(c) for _, c in changed)}")

    if conflicts:
        print()
        print("=== 手填值被覆盖的基因（需审阅）===")
        for name, gene, old, new in conflicts:
            print(f"  {name}: {gene} {old} -> {new}")

    if dry_run:
        print("\n[dry-run] 未写回文件，加 --dry-run 预览")


if __name__ == "__main__":
    main()
