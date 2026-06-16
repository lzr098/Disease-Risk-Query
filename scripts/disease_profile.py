"""DiseaseProfile data model.

This module defines the central abstraction for disease risk analysis:
a structured, serializable description of a disease's genomic analysis
space (genes, regulatory regions, known variants) and its contribution
model (how variants are weighted toward disease risk).
"""

from __future__ import annotations

import json
import logging
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Optional

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class GeneWeight:
    """A gene linked to the disease with evidence and weight."""

    gene: str
    tier: str  # mendelian_high | mendelian_mod | strong_gwas | gwas
    contribution_score: float
    penetrance: str  # >0.95 | high | moderate | low_to_moderate | low | dosage | very_low
    penetrance_score: float
    evidence: str  # familial | rare_variant | gwas | mixed | common_risk
    note: str = ""
    phenotype_assoc: str = ""  # Gene-disease phenotype association from MyGene/OMIM
    key_domains: str = ""  # Key protein domains from MyGene InterPro / UniProt
    clingen_validity: str = ""  # ClinGen validity: Definitive / Strong / Moderate / Limited / Disputed / Refuted
    is_mendelian: bool = False  # True if tier is mendelian_high or mendelian_mod

    def __post_init__(self):
        # Auto-derive is_mendelian from tier
        if self.tier in ("mendelian_high", "mendelian_mod"):
            object.__setattr__(self, "is_mendelian", True)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "GeneWeight":
        filtered = {k: v for k, v in data.items() if k in cls.__dataclass_fields__}
        # Auto-derive is_mendelian from tier if not explicitly set
        if "is_mendelian" not in filtered or not filtered["is_mendelian"]:
            tier = filtered.get("tier", "")
            filtered["is_mendelian"] = tier in ("mendelian_high", "mendelian_mod")
        return cls(**filtered)


@dataclass(frozen=True)
class VariantWeight:
    """A known variant (pathogenic, GWAS lead, PRS, dosage risk)."""

    chrom: str
    pos: int
    ref: Optional[str] = None
    alt: Optional[str] = None
    rsid: Optional[str] = None
    gene: Optional[str] = None
    effect_allele: Optional[str] = None
    other_allele: Optional[str] = None
    beta: Optional[float] = None
    or_value: Optional[float] = None
    eaf_eur: Optional[float] = None
    eaf_eas: Optional[float] = None
    variant_class: str = "gwas_lead"  # known_pathogenic | gwas_lead | prs | dosage_risk
    contribution_score: float = 0.0
    confidence: str = "moderate"
    note: str = ""

    def __post_init__(self):
        # Infer ref/alt from effect/other alleles if missing
        if not self.ref or not self.alt:
            object.__setattr__(
                self, "ref",
                (self.ref or self.other_allele or "").upper(),
            )
            object.__setattr__(
                self, "alt",
                (self.alt or self.effect_allele or "").upper(),
            )

    @property
    def vcf_key(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}:{self.alt}"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "VariantWeight":
        return cls(**{k: v for k, v in data.items() if k in cls.__dataclass_fields__})


@dataclass(frozen=True)
class RegulatoryRegion:
    """A regulatory region linked to a disease gene."""

    chrom: str
    start: int  # 0-based
    end: int
    gene: Optional[str] = None
    rtype: str = "regulatory"  # promoter | enhancer | utr | regulatory
    source: str = "literature"
    note: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "RegulatoryRegion":
        return cls(**{k: v for k, v in data.items() if k in cls.__dataclass_fields__})


@dataclass(frozen=True)
class DomainRegion:
    name: str
    residues: str
    note: str = ""

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "DomainRegion":
        return cls(**{k: v for k, v in data.items() if k in cls.__dataclass_fields__})


@dataclass(frozen=True)
class CriticalResidue:
    residue: str
    note: str = ""

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "CriticalResidue":
        return cls(**{k: v for k, v in data.items() if k in cls.__dataclass_fields__})


@dataclass(frozen=True)
class GeneDomainInfo:
    note: str = ""
    regions: list[DomainRegion] = field(default_factory=tuple)
    critical_residues: list[CriticalResidue] = field(default_factory=tuple)

    def __post_init__(self):
        # Frozen dataclass: convert lists to tuples for hashability if needed.
        object.__setattr__(self, "regions", tuple(self.regions))
        object.__setattr__(self, "critical_residues", tuple(self.critical_residues))

    def to_dict(self) -> dict[str, Any]:
        return {
            "note": self.note,
            "regions": [asdict(r) for r in self.regions],
            "critical_residues": [asdict(r) for r in self.critical_residues],
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "GeneDomainInfo":
        return cls(
            note=data.get("note", ""),
            regions=[DomainRegion.from_dict(r) for r in data.get("regions", [])],
            critical_residues=[CriticalResidue.from_dict(r) for r in data.get("critical_residues", [])],
        )


@dataclass
class DiseaseProfile:
    """Central description of a disease for risk analysis."""

    canonical_name: str
    aliases: list[str] = field(default_factory=list)
    mode: str = "mendelian"  # mendelian | complex | mixed
    hpo_terms: list[dict] = field(default_factory=list)
    gene_set: list[GeneWeight] = field(default_factory=list)
    known_pathogenic_variants: list[VariantWeight] = field(default_factory=list)
    gwas_lead_snps: list[VariantWeight] = field(default_factory=list)
    prs_variants: list[VariantWeight] = field(default_factory=list)
    prs_variants_high: list[VariantWeight] = field(default_factory=list)
    regulatory_regions: list[RegulatoryRegion] = field(default_factory=list)
    key_regions: dict[str, GeneDomainInfo] = field(default_factory=dict)
    key_literature: list[dict] = field(default_factory=list)
    contribution_model: dict = field(default_factory=dict)

    @property
    def all_genes(self) -> set[str]:
        return {g.gene for g in self.gene_set}

    @property
    def all_known_variants(self) -> list[VariantWeight]:
        return (
            self.known_pathogenic_variants
            + self.gwas_lead_snps
            + self.prs_variants
            + self.prs_variants_high
        )

    @property
    def core_genes(self) -> set[str]:
        return {g.gene for g in self.gene_set if g.tier in ("mendelian_high", "mendelian_mod")}

    @property
    def gwas_loci_genes(self) -> set[str]:
        return {g.gene for g in self.gene_set if g.tier not in ("mendelian_high", "mendelian_mod")}

    @property
    def gene_contribution_map(self) -> dict[str, float]:
        return {g.gene: g.contribution_score for g in self.gene_set}

    @property
    def gene_penetrance_map(self) -> dict[str, str]:
        return {g.gene: g.penetrance for g in self.gene_set}

    @property
    def gene_tier_map(self) -> dict[str, str]:
        return {g.gene: g.tier for g in self.gene_set}

    @property
    def gene_phenotype_map(self) -> dict[str, str]:
        """Gene → phenotype_assoc description."""
        return {g.gene: g.phenotype_assoc for g in self.gene_set if g.phenotype_assoc}

    @property
    def gene_domain_map(self) -> dict[str, str]:
        """Gene → key_domains annotations."""
        return {g.gene: g.key_domains for g in self.gene_set if g.key_domains}

    @property
    def gene_clingen_map(self) -> dict[str, str]:
        """Gene → clingen_validity rating."""
        return {g.gene: g.clingen_validity for g in self.gene_set if g.clingen_validity}

    @property
    def snp_contribution_map(self) -> dict[str, float]:
        result: dict[str, float] = {}
        for v in self.all_known_variants:
            key = v.rsid or v.vcf_key
            result[key] = v.contribution_score
        return result

    def to_dict(self) -> dict[str, Any]:
        return {
            "canonical_name": self.canonical_name,
            "aliases": self.aliases,
            "mode": self.mode,
            "hpo_terms": self.hpo_terms,
            "gene_set": [g.to_dict() for g in self.gene_set],
            "known_pathogenic_variants": [v.to_dict() for v in self.known_pathogenic_variants],
            "gwas_lead_snps": [v.to_dict() for v in self.gwas_lead_snps],
            "prs_variants": [v.to_dict() for v in self.prs_variants],
            "prs_variants_high": [v.to_dict() for v in self.prs_variants_high],
            "regulatory_regions": [r.to_dict() for r in self.regulatory_regions],
            "key_regions": {k: v.to_dict() for k, v in self.key_regions.items()},
            "key_literature": self.key_literature,
            "contribution_model": self.contribution_model,
        }

    def to_json(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2, ensure_ascii=False))

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "DiseaseProfile":
        return cls(
            canonical_name=data["canonical_name"],
            aliases=list(data.get("aliases", [])),
            mode=data.get("mode", "mendelian"),
            hpo_terms=list(data.get("hpo_terms", [])),
            gene_set=[GeneWeight.from_dict(g) for g in data.get("gene_set", [])],
            known_pathogenic_variants=[VariantWeight.from_dict(v) for v in data.get("known_pathogenic_variants", [])],
            gwas_lead_snps=[VariantWeight.from_dict(v) for v in data.get("gwas_lead_snps", [])],
            prs_variants=[VariantWeight.from_dict(v) for v in data.get("prs_variants", [])],
            prs_variants_high=[VariantWeight.from_dict(v) for v in data.get("prs_variants_high", [])],
            regulatory_regions=[RegulatoryRegion.from_dict(r) for r in data.get("regulatory_regions", [])],
            key_regions={k: GeneDomainInfo.from_dict(v) for k, v in data.get("key_regions", {}).items()},
            key_literature=list(data.get("key_literature", [])),
            contribution_model=dict(data.get("contribution_model", {})),
        )

    @classmethod
    def from_json(cls, path: Path) -> "DiseaseProfile":
        data = json.loads(path.read_text(encoding="utf-8"))
        return cls.from_dict(data)


def build_default_contribution_model(mode: str) -> dict[str, Any]:
    """Return a default contribution model for the disease mode."""
    base = {
        "mendelian_high": {
            "scoring": "binary_presence",
            "weight": 1.0,
            "description": "High-penetrance pathogenic variant; presence implies high risk.",
        },
        "mendelian_mod": {
            "scoring": "gene_contribution_x_penetrance_x_zygosity",
            "weight": 0.8,
            "description": "Moderate-penetrance variant in disease gene.",
        },
        "known_pathogenic": {
            "scoring": "template_curated_pathogenic_x_zygosity_x_penetrance",
            "weight": 0.9,
            "description": (
                "Curated known pathogenic variant from disease template; "
                "scored by dosage, gene penetrance, and template contribution score."
            ),
        },
        "dosage_risk": {
            "scoring": "copy_number_x_log_or",
            "weight": 0.5,
            "description": "Risk allele with dosage effect (e.g. APOE e4).",
        },
        "gwas_prs": {
            "scoring": "weighted_sum_sqrt_beta_x_dosage",
            "weight": 0.3,
            "description": (
                "Common variant polygenic contribution; "
                "uses sqrt(|beta|) to amplify small-effect variants."
            ),
        },
        "prs_high": {
            "scoring": "high_confidence_prs_weighted_sum",
            "weight": 0.9,
            "description": (
                "High-confidence PRS variants from validated "
                "polygenic risk scores with published weights."
            ),
        },
        "regulatory": {
            "scoring": "gene_contribution_x_predicted_impact",
            "weight": 0.1,
            "description": "Variant in regulatory region of core gene (weak evidence).",
        },
    }
    if mode == "complex":
        # In complex mode, GWAS/PRS layer is the primary signal.
        base["gwas_prs"]["weight"] = 0.8
        base["mendelian_high"]["weight"] = 0.7
    return base
