"""Tests for DiseaseProfile-driven refactor modules."""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

# Ensure skill scripts are importable
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

from disease_profile import (
    DiseaseProfile,
    GeneDomainInfo,
    GeneWeight,
    RegulatoryRegion,
    VariantWeight,
    build_default_contribution_model,
)
from disease_profile_builder import _normalize_builtin_ref, build_or_load_profile
from disease_space_query import _build_unified_bed, query_known_variants
from contribution_scorer import score
from constants import DISEASE_BUILTIN_REFS
from clinvar_phenotype_matcher import (
    classify_clinvar_sig,
    enrich_variants_with_clinvar,
    is_clinvar_plp,
    is_clinvar_vus,
)
from literature_source import build_literature_for_profile


def test_variant_weight_vcf_key_with_ref_alt():
    v = VariantWeight(chrom="chr19", pos=44908822, ref="C", alt="T", rsid="rs429358")
    assert v.vcf_key == "chr19:44908822:C:T"


def test_variant_weight_infers_ref_alt_from_effect_alleles():
    v = VariantWeight(chrom="chr19", pos=44908822, effect_allele="T", other_allele="C")
    assert v.ref == "C"
    assert v.alt == "T"


def test_disease_profile_roundtrip():
    profile = DiseaseProfile(
        canonical_name="test disease",
        aliases=["test"],
        mode="mendelian",
        gene_set=[
            GeneWeight(
                gene="GENE1", tier="mendelian_high",
                contribution_score=1.0, penetrance=">0.95",
                penetrance_score=0.95, evidence="familial",
            ),
        ],
        gwas_lead_snps=[
            VariantWeight(
                chrom="chr1", pos=1000, ref="A", alt="G",
                rsid="rs123", effect_allele="G", beta=0.1,
            ),
        ],
        regulatory_regions=[
            RegulatoryRegion(chrom="chr1", start=900, end=1100, gene="GENE1"),
        ],
        key_regions={
            "GENE1": GeneDomainInfo(
                note="test domain",
                regions=[],
                critical_residues=[],
            ),
        },
        contribution_model=build_default_contribution_model("mendelian"),
    )
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "profile.json"
        profile.to_json(path)
        loaded = DiseaseProfile.from_json(path)
    assert loaded.canonical_name == "test disease"
    assert loaded.all_genes == {"GENE1"}
    assert len(loaded.all_known_variants) == 1
    assert len(loaded.regulatory_regions) == 1
    assert "GENE1" in loaded.key_regions


def test_normalize_builtin_ref_includes_all_required_fields():
    for key in DISEASE_BUILTIN_REFS:
        normalized = _normalize_builtin_ref(key, DISEASE_BUILTIN_REFS[key])
        for field in [
            "canonical_name", "aliases", "mode", "hpo_terms", "gene_set",
            "known_pathogenic_variants", "gwas_lead_snps", "prs_variants",
            "regulatory_regions", "key_regions", "key_literature", "contribution_model",
        ]:
            assert field in normalized, f"{key} missing {field}"
        assert normalized["contribution_model"]
        assert isinstance(normalized["gene_set"], list)
        assert isinstance(normalized["regulatory_regions"], list)


def test_builtin_profiles_have_consistent_runtime_schema():
    for disease in ["alzheimer disease", "parkinson disease", "adult vision disorders", "hyperuricemia"]:
        profile = build_or_load_profile(disease, refresh=True)
        assert profile.canonical_name
        assert profile.mode in ("mendelian", "complex", "mixed")
        assert profile.gene_set
        # All profiles should have a contribution model
        assert profile.contribution_model
        # key_regions should be dict[str, dict]
        assert isinstance(profile.key_regions, dict)


def test_build_unified_bed():
    profile = DiseaseProfile(
        canonical_name="test",
        gene_set=[
            GeneWeight(
                gene="GENE1", tier="mendelian_high",
                contribution_score=1.0, penetrance=">0.95",
                penetrance_score=0.95, evidence="familial",
            ),
        ],
        regulatory_regions=[
            RegulatoryRegion(chrom="chr1", start=500, end=600, gene="GENE1"),
        ],
    )
    gene_coords = {
        "GENE1": [("chr1", 100, 200, "+")],
    }
    with tempfile.TemporaryDirectory() as tmp:
        bed_path = Path(tmp) / "unified.bed"
        count = _build_unified_bed(profile, gene_coords, bed_path)
        assert count == 2
        lines = bed_path.read_text().strip().split("\n")
        assert lines[0].startswith("chr1\t100\t200\tGENE1")
        assert lines[1].startswith("chr1\t500\t600\tGENE1")


def test_query_known_variants_with_ref_ref_inference():
    header = [
        "##fileformat=VCFv4.2",
        "##reference=GRCh38",
        "##contig=<ID=chr21,length=46709983>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
    ]
    records = [
        "chr21\t25891784\t.\tC\tA\t99\tPASS\t.\tGT\t0/1",
    ]
    variants = [
        VariantWeight(chrom="chr21", pos=25891784, ref="C", alt="A", effect_allele="A"),
        VariantWeight(chrom="chr21", pos=25891785, ref="G", alt="T", effect_allele="T"),
    ]
    with tempfile.TemporaryDirectory() as tmp:
        src = Path(tmp) / "in.vcf.gz"
        _write_vcf(src, header, records)
        result = query_known_variants(src, variants)
    assert len(result) == 2
    found = [r for r in result if not r.inferred_ref_ref]
    inferred = [r for r in result if r.inferred_ref_ref]
    assert len(found) == 1
    assert len(inferred) == 1
    assert found[0].gt == "0/1"
    assert found[0].dosage == 1
    assert inferred[0].gt == "0/0"
    assert inferred[0].dosage == 0


def test_contribution_scorer_mendelian_high_dominates():
    profile = DiseaseProfile(
        canonical_name="test",
        mode="mendelian",
        gene_set=[
            GeneWeight(
                gene="PSEN1", tier="mendelian_high",
                contribution_score=1.0, penetrance=">0.95",
                penetrance_score=0.98, evidence="familial",
            ),
        ],
        contribution_model=build_default_contribution_model("mendelian"),
    )
    tiered = {
        "tier1_variants": [
            {"gene": "PSEN1", "chrom": "chr14", "pos": 73173707, "gt": "0/1", "impact": "HIGH"},
        ],
        "tier2_variants": [],
        "tier3_variants": [],
    }
    result = score(profile, tiered, [], {})
    assert result.mendelian_high
    assert result.overall_level == "high"
    assert result.overall_score >= 0.8


def test_contribution_scorer_tier3_pathogenic_core_gene_scored():
    """Tier 3 pathogenic variants in core genes must be scored (not ignored)."""
    profile = DiseaseProfile(
        canonical_name="test",
        mode="mendelian",
        gene_set=[
            GeneWeight(
                gene="APP", tier="mendelian_high",
                contribution_score=1.0, penetrance=">0.95",
                penetrance_score=0.98, evidence="familial",
            ),
        ],
        contribution_model=build_default_contribution_model("mendelian"),
    )
    tiered = {
        "tier1_variants": [],
        "tier2_variants": [],
        "tier3_variants": [
            {
                "gene": "APP", "chrom": "chr21", "pos": 25891784,
                "ref": "C", "alt": "A", "gt": "0/1",
                "impact": "HIGH", "clinvar_sig": "Pathogenic",
            },
        ],
    }
    result = score(profile, tiered, [], {})
    assert result.mendelian_mod
    app_hit = next(h for h in result.mendelian_mod if h["gene"] == "APP")
    assert app_hit["tier"] == 3
    assert app_hit["contribution"] > 0


def test_contribution_scorer_dosage_risk():
    profile = DiseaseProfile(
        canonical_name="test",
        mode="mendelian",
        gwas_lead_snps=[
            VariantWeight(
                chrom="chr19", pos=44908822, effect_allele="T", other_allele="C",
                or_value=3.0, variant_class="dosage_risk", contribution_score=0.5,
            ),
        ],
        contribution_model=build_default_contribution_model("mendelian"),
    )
    from disease_space_query import KnownVariantGenotype
    kg = KnownVariantGenotype(
        variant=profile.gwas_lead_snps[0],
        chrom="chr19", pos=44908822, ref="C", alt="T",
        gt="0/1", dosage=1, inferred_ref_ref=False,
        filter_status="PASS", sample_format="GT", sample_values="0/1",
    )
    result = score(profile, {"tier1_variants": [], "tier2_variants": [], "tier3_variants": []}, [kg], {})
    assert len(result.dosage_risk) == 1
    assert result.dosage_risk[0]["dosage"] == 1


def _write_vcf(path: Path, header: list[str], records: list[str]) -> None:
    import subprocess
    raw = "\n".join(header + records) + "\n"
    if str(path).endswith(".gz"):
        proc = subprocess.run(
            ["bcftools", "view", "-Oz", "-o", str(path)],
            input=raw, text=True, check=True,
        )
        subprocess.run(["bcftools", "index", str(path)], check=True)
    else:
        path.write_text(raw, encoding="utf-8")


def _write_clinvar_vcf(path: Path, records: list[str]) -> None:
    """Write a minimal ClinVar-style VCF for testing."""
    import subprocess
    header = [
        "##fileformat=VCFv4.1",
        "##reference=GRCh38",
        "##contig=<ID=chr1,length=248956422>",
        '##INFO=<ID=CLNDN,Number=.,Type=String,Description="ClinVar disease name">',
        '##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance">',
        '##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="Review status">',
        '##INFO=<ID=ALLELEID,Number=1,Type=Integer,Description="Allele ID">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    raw = "\n".join(header + records) + "\n"
    proc = subprocess.run(
        ["bcftools", "view", "-Oz", "-o", str(path)],
        input=raw, text=True, check=True,
    )
    subprocess.run(["bcftools", "index", str(path)], check=True)


def test_classify_clinvar_sig():
    assert classify_clinvar_sig("Pathogenic") == "pathogenic"
    assert classify_clinvar_sig("Likely_pathogenic") == "likely_pathogenic"
    assert classify_clinvar_sig("Uncertain_significance") == "vus"
    assert classify_clinvar_sig("Conflicting_interpretations_of_pathogenicity") == "conflicting"
    assert classify_clinvar_sig("risk_factor") == "risk_factor"
    assert classify_clinvar_sig("protective") == "protective"
    assert is_clinvar_plp("Pathogenic")
    assert is_clinvar_plp("Likely_pathogenic")
    assert not is_clinvar_plp("Uncertain_significance")
    assert is_clinvar_vus("Uncertain_significance")


def test_enrich_variants_with_clinvar_does_not_filter_vus():
    """VUS variants must be enriched but never dropped."""
    records = [
        "chr1\t1000\t.\tA\tG\t.\t.\tCLNDN=hyperuricemia|not_provided;CLNSIG=Uncertain_significance;CLNREVSTAT=criteria_provided;ALLELEID=1",
        "chr1\t2000\t.\tC\tT\t.\t.\tCLNDN=hyperuricemia;CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided;ALLELEID=2",
    ]
    variants = [
        {"chrom": "chr1", "pos": 1000, "ref": "A", "alt": "G"},
        {"chrom": "chr1", "pos": 2000, "ref": "C", "alt": "T"},
    ]
    with tempfile.TemporaryDirectory() as tmp:
        clinvar_path = Path(tmp) / "clinvar.vcf.gz"
        _write_clinvar_vcf(clinvar_path, records)
        enriched = enrich_variants_with_clinvar(variants, "hyperuricemia", vcf_path=clinvar_path)

    assert len(enriched) == 2
    vus = next(v for v in enriched if v["pos"] == 1000)
    plp = next(v for v in enriched if v["pos"] == 2000)
    assert vus["clinvar_category"] == "vus"
    assert vus["clinvar_plp"] is False
    assert plp["clinvar_category"] == "pathogenic"
    assert plp["clinvar_plp"] is True
    assert plp["clinvar_disease_match"] is True


def test_enrich_variants_with_clinvar_conflicting_preserved():
    records = [
        "chr1\t3000\t.\tG\tA\t.\t.\tCLNDN=hyperuricemia;CLNSIG=Conflicting_interpretations_of_pathogenicity;CLNREVSTAT=criteria_provided;ALLELEID=3",
    ]
    variants = [{"chrom": "chr1", "pos": 3000, "ref": "G", "alt": "A"}]
    with tempfile.TemporaryDirectory() as tmp:
        clinvar_path = Path(tmp) / "clinvar.vcf.gz"
        _write_clinvar_vcf(clinvar_path, records)
        enriched = enrich_variants_with_clinvar(variants, "hyperuricemia", vcf_path=clinvar_path)
    assert len(enriched) == 1
    assert enriched[0]["clinvar_category"] == "conflicting"
    assert enriched[0]["clinvar_plp"] is False


def test_build_literature_for_profile_combines_curated_and_dynamic(monkeypatch):
    """Curated entries are preserved; dynamic PubMed entries are appended."""
    monkeypatch.setattr(
        "literature_source._pubmed_query",
        lambda disease, genes, retmax: [
            {
                "pmid": "12345",
                "title": "Test PubMed article",
                "authors": ["Smith J"],
                "journal": "Test Journal",
                "year": "2024",
                "genes": [],
                "note": "PubMed search",
                "source": "pubmed",
            },
        ],
    )
    with tempfile.TemporaryDirectory() as tmp:
        lit = build_literature_for_profile(
            "hyperuricemia",
            ["SLC2A9", "ABCG2"],
            cache_root=Path(tmp),
            use_pubmed=True,
        )
    assert len(lit) >= 1
    sources = {e.get("source") for e in lit}
    assert "pubmed" in sources


def test_build_literature_cache_hit(monkeypatch):
    """Second call with same cache should not re-query PubMed."""
    call_count = {"n": 0}

    def fake_query(disease, genes, retmax):
        call_count["n"] += 1
        return [{"pmid": "99999", "title": "Cached", "source": "pubmed", "genes": []}]

    monkeypatch.setattr("literature_source._pubmed_query", fake_query)
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        lit1 = build_literature_for_profile("test disease", ["GENE1"], cache_root=root)
        lit2 = build_literature_for_profile("test disease", ["GENE1"], cache_root=root)
    assert lit1 == lit2
    assert call_count["n"] == 1
