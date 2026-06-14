"""Basic unit tests for disease risk query modules."""

from __future__ import annotations

import gzip
import sys
import tempfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

from hpo_mapper import resolve_disease_query
from gene_set_builder import build_disease_gene_set
from risk_scorer import calculate_total_score
from clinvar_phenotype_matcher import _normalize_text


def _write_vcf(path: Path, header_lines: list[str], records: list[str]) -> None:
    with gzip.open(path, "wt") as f:
        for line in header_lines:
            f.write(line + "\n")
        for line in records:
            f.write(line + "\n")



def test_hpo_mapper_alzheimer():
    result = resolve_disease_query("Alzheimer disease")
    assert result["hpo_id"] == "HP:0000726"
    assert result["matched_by"] in ("curated", "curated_builtin")
    assert "APOE" in result["genes"] or "APP" in result["genes"]


def test_hpo_mapper_hyperuricemia_chinese():
    result = resolve_disease_query("高尿酸血症")
    assert result["hpo_id"] == "HP:0002149"
    assert result["hpo_name"] == "Hyperuricemia"
    assert result["matched_by"] == "curated_builtin"


def test_gene_set_builder():
    hpo = resolve_disease_query("Alzheimer disease")
    gs = build_disease_gene_set(
        "Alzheimer disease",
        hpo_genes=hpo["genes"],
        max_genes=20,
    )
    assert gs["total"] <= 20
    assert "APOE" in gs["merged_genes"] or "APP" in gs["merged_genes"]
    # No non-gene symbols should leak through
    for gene in gs["merged_genes"]:
        assert not gene.isdigit()
        assert len(gene) >= 2


def test_risk_scorer_moderate_with_tier1():
    tier1 = [
        {"CHROM": "21", "POS": 100, "REF": "A", "ALT": "G", "GENE": "APP", "GT": "1/1", "IMPACT": "HIGH", "Consequence": "stop_gained", "gnomAD_AF": 0.0001},
        {"CHROM": "14", "POS": 200, "REF": "C", "ALT": "T", "GENE": "PSEN1", "GT": "0/1", "IMPACT": "HIGH", "Consequence": "frameshift_variant", "gnomAD_AF": 0.0001},
    ]
    score = calculate_total_score(tier1, [], [], [], set(), 5)
    assert 50 <= score["total_score"] <= 79
    assert score["risk_level"] == "中等风险"


def test_risk_scorer_high():
    tier1 = [
        {"CHROM": "21", "POS": 100, "REF": "A", "ALT": "G", "GENE": "APP", "GT": "1/1", "IMPACT": "HIGH", "Consequence": "stop_gained", "gnomAD_AF": 0.0001},
    ]
    tier2 = [
        {"CHROM": "1", "POS": 300, "REF": "G", "ALT": "A", "GENE": "PSEN2", "GT": "0/1", "IMPACT": "MODERATE", "Consequence": "missense_variant", "gnomAD_AF": 0.0001},
    ]
    lit_variants = [{"CHROM": "21", "POS": 100, "REF": "A", "ALT": "G"}]
    score = calculate_total_score(tier1, tier2, ["PSEN2"], lit_variants, {"PSEN2"}, 5)
    assert score["total_score"] >= 80
    assert score["risk_level"] == "高风险"


def test_risk_scorer_empty():
    score = calculate_total_score([], [], [], [], set(), 0)
    assert score["total_score"] < 20
    assert score["risk_level"] == "无明确风险"


def test_liftover_grch37_to_grch38():
    from liftover import liftover_vcf

    header = [
        "##fileformat=VCFv4.2",
        "##reference=GRCh37",
        "##contig=<ID=chr14,length=107349540>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
    ]
    records = [
        "chr14\t73145000\t.\tC\tT\t99\tPASS\t.\tGT\t0/1",
        "chr14\t73152000\t.\tA\tG\t99\tPASS\t.\tGT\t1/1",
    ]
    with tempfile.TemporaryDirectory() as tmp:
        src = Path(tmp) / "in.vcf.gz"
        dst = Path(tmp) / "out.vcf.gz"
        _write_vcf(src, header, records)

        result = liftover_vcf(src, dst, source_build="GRCh37", target_build="GRCh38")
        assert result["lifted"] is True
        assert result["input_records"] == 2
        assert result["output_records"] >= 1
        assert Path(result["output_path"]).exists()
        assert result["validation"]["pass"] is True


def test_vep115_compat_annotation_and_parser():
    from gpa_compat import patch_gpa_csq_parser, pre_annotate_with_vep

    header = [
        "##fileformat=VCFv4.2",
        "##reference=GRCh38",
        "##contig=<ID=chr21,length=46709983>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
    ]
    records = ["chr21\t25891784\t.\tC\tA\t99\tPASS\t.\tGT\t0/1"]
    with tempfile.TemporaryDirectory() as tmp:
        src = Path(tmp) / "in.vcf.gz"
        dst = Path(tmp) / "vep115.vcf.gz"
        _write_vcf(src, header, records)

        annotated = pre_annotate_with_vep(src, dst)
        assert Path(annotated).exists()

        # Header should contain modern gnomAD field names
        csq = None
        with gzip.open(annotated, "rt") as f:
            for line in f:
                if line.startswith("##INFO=<ID=CSQ"):
                    csq = line
                    break
        assert csq is not None
        assert "gnomADe_AF" in csq or "gnomADg_AF" in csq

        # Patched parser should populate gnomAD_AF
        sys.path.insert(0, str(Path.home() / ".workbuddy" / "skills" / "dgra-genomic-risk" / "scripts"))
        patch_gpa_csq_parser()
        from gpa_input import parse_annotated_vcf
        variants = parse_annotated_vcf(str(annotated))
        assert len(variants) == 1
        assert variants[0]["GENE"] == "APP"
        assert variants[0]["IMPACT"] in {"MODERATE", "HIGH"}


def test_ad_template_excludes_non_ad_genes():
    from disease_reference import get_disease_reference
    ref = get_disease_reference("alzheimer disease")
    gene_set = {g["gene"] for g in ref["gene_set"]}
    non_ad = {"VWF", "CD36", "MYH11"}
    assert non_ad.isdisjoint(gene_set), f"Non-AD genes leaked into AD template: {non_ad & gene_set}"


def test_clinvar_disease_filter_rejects_non_ad_pathogenic():
    from clinvar_phenotype_matcher import filter_variants_by_clinvar_disease

    variants = [
        {"CHROM": "chr12", "POS": 6019487, "REF": "G", "ALT": "A", "GENE": "VWF", "clinvar_diseases": ["Hereditary_von_Willebrand_disease"]},
        {"CHROM": "chr7", "POS": 80661109, "REF": "AAC", "ALT": "A", "GENE": "CD36", "clinvar_diseases": ["CD36-related_disorder"]},
        {"CHROM": "chr16", "POS": 15735448, "REF": "G", "ALT": "A", "GENE": "MYH11", "clinvar_diseases": ["Aortic_aneurysm"]},
    ]
    filtered = filter_variants_by_clinvar_disease(variants, "alzheimer disease", require_match=True)
    assert len(filtered) == 0, "Non-AD pathogenic variants must be filtered out when require_match=True"


def test_omim_alzheimer_english_keywords_and_clean_symbols():
    from gene_set_builder import build_disease_gene_set

    # Chinese disease name but English OMIM keyword (as pipeline now supplies)
    result = build_disease_gene_set(
        disease_name="阿尔茨海默病",
        omim_keywords=["alzheimer disease"],
        max_genes=200,
    )
    omim_genes = set(result["omim_genes"])
    # Core AD genes should be present
    assert "APP" in omim_genes
    assert "PSEN1" in omim_genes
    # Disease names and locus labels must be cleaned out
    assert "ALZHEIMER DISEASE" not in omim_genes
    assert not any("\n" in g for g in omim_genes), "OMIM symbols should not contain embedded newlines"


def test_clinvar_disease_keywords_chinese_no_empty():
    from clinvar_phenotype_matcher import disease_keywords

    kws = disease_keywords("高尿酸血症")
    assert "hyperuricemia" in kws
    assert "gout" in kws
    assert "" not in kws
    assert not any(_normalize_text(kw) == "" for kw in kws)


def test_resolve_disease_mode_unknown_disease():
    from constants import resolve_disease_mode, DISEASE_MODE_COMPLEX, DISEASE_MODE_MENDELIAN

    # Unknown complex-trait keyword should resolve to complex mode
    assert resolve_disease_mode("diabetes") == DISEASE_MODE_COMPLEX
    assert resolve_disease_mode("hypertension") == DISEASE_MODE_COMPLEX
    # Unknown disease without complex keyword defaults to mendelian
    assert resolve_disease_mode("foobar syndrome") == DISEASE_MODE_MENDELIAN


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
