"""Constants and default paths for gpa-disease-risk-query skill."""

from __future__ import annotations

from pathlib import Path

# Base directories
WORKBUDDY_DATA = Path.home() / ".workbuddy" / "data"
WORKBUDDY_SKILL = Path(__file__).resolve().parent.parent

# Local bioinformatics assets (reused from other skills)
GRCH38_FASTA = WORKBUDDY_DATA / "genome" / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
CLINVAR_VCF = WORKBUDDY_DATA / "clinvar" / "clinvar.vcf.gz"
OMIM_DB = WORKBUDDY_DATA / "omim" / "omim.db"
HPO_GENES_TO_PHENOTYPE = WORKBUDDY_DATA / "hpo" / "genes_to_phenotype.txt"
HGNC_LOOKUP = WORKBUDDY_DATA / "hgnc" / "hgnc_lookup.json"
GENCODE_GTF = WORKBUDDY_DATA / "gencode" / "gencode.v44.annotation.gtf.gz"

# dgra-prefilter reference BEDs (used for gene coordinate extraction)
DGRA_PREFILTER_REFS = Path.home() / ".dgra-prefilter" / "refs"
GENCODE_GENE_LOCI_BED = DGRA_PREFILTER_REFS / "gencode_v44_gene_loci.bed"
CLINVAR_PATHOGENIC_BED = DGRA_PREFILTER_REFS / "clinvar_pathogenic_GRCh38.bed"

# LiftOver cache (UCSC chain files via pyliftover)
LIFTOVER_DIR = WORKBUDDY_DATA / "liftover"
HG19_TO_HG38_CHAIN = LIFTOVER_DIR / "hg19ToHg38.over.chain.gz"
HG38_TO_HG19_CHAIN = LIFTOVER_DIR / "hg38ToHg19.over.chain.gz"
LIFTOVER_URLS = {
    ("hg19", "hg38"): "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
    ("hg38", "hg19"): "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
}

# VEP Docker defaults
VEP_IMAGE = "ensemblorg/ensembl-vep:latest"
VEP_CACHE = Path.home() / ".workbuddy" / "tools" / "vep" / "cache"

# Skill-internal cache / references
SKILL_REFS = WORKBUDDY_SKILL / "references"
GENE_COORDS_CACHE = SKILL_REFS / "gencode_v44_gene_coords.json"
HPO_DISEASE_CACHE = SKILL_REFS / "hpo_disease_terms.json"

# Derived cache from OMIM DB to speed up disease-gene mapping
OMIM_TITLE_GENE_CACHE = WORKBUDDY_DATA / "omim" / "omim_title_gene_cache.json"

# Default outputs
DEFAULT_OUTPUT_DIR = Path.cwd()

# Supported genome builds
SUPPORTED_BUILDS = {"GRCh38", "GRCh37", "hg19", "hg38"}
GRCH38_ALIASES = {"GRCh38", "hg38"}
GRCH37_ALIASES = {"GRCh37", "hg19"}

# Tissue mapping for common disease categories
DEFAULT_TISSUE = "general"
DISEASE_TISSUE_MAP = {
    "neurological": [
        "alzheimer", "alzheimer disease", "parkinson", "epilepsy", "seizure",
        "neurodegenerative", "dementia", "amyotrophic lateral sclerosis", "als",
        "huntington", "ataxia", "neuropathy", "migraine", "autism", "schizophrenia",
        "bipolar", "intellectual disability", "developmental delay", "encephalopathy",
    ],
    "cardiovascular": [
        "cardiomyopathy", "arrhythmia", "long qt", "brugada", "hypertrophic cardiomyopathy",
        "dilated cardiomyopathy", "congenital heart", "aortic", "marfan", "loeydietz",
        "vascular", "stroke", "coronary", "myocardial",
    ],
    "renal": [
        "polycystic kidney", "nephrotic syndrome", "nephritis", "renal", "kidney",
        "alport", "fabry", "nephronophthisis",
    ],
    "hepatic": [
        "liver", "hepatic", "cirrhosis", "hemochromatosis", "wilson", "alpha-1 antitrypsin",
    ],
    "hematopoietic": [
        "anemia", "thalassemia", "sickle cell", "hemophilia", "leukemia", "lymphoma",
        "myelodysplastic", "bone marrow", "bleeding disorder", "coagulation",
    ],
}

# Risk score weights (must match PRD section 6A)
WEIGHT_TIER1 = 40
WEIGHT_TIER2 = 25
WEIGHT_LITERATURE = 20
WEIGHT_RARITY = 10
WEIGHT_SEX_AGE = 5

# Tier contribution sub-rules
TIER1_BASE = 40
TIER1_HOM_PURE_LOF_BONUS = 20
TIER2_BASE = 15
TIER2_MULTIHIT_BONUS = 10
LIT_VARIANT_BONUS = 15
LIT_GENE_BONUS = 5

# Risk thresholds (mendelian mode)
RISK_HIGH_THRESHOLD = 80
RISK_MODERATE_THRESHOLD = 50
RISK_LOW_THRESHOLD = 20

# Complex-trait contribution score weights (0-100)
# Used when disease mode is "complex" instead of the tier-based risk score.
COMPLEX_WEIGHT_MONOGENIC = 30          # high-penetrance single-gene variants (Tier 1/2)
COMPLEX_WEIGHT_RARE_FUNCTIONAL = 25    # rare missense/LoF in core genes (Tier 3 filtered)
COMPLEX_WEIGHT_GWAS_COMMON = 25        # direct GWAS lead SNP genotypes
COMPLEX_WEIGHT_LITERATURE = 10         # literature-supported genes/variants
COMPLEX_WEIGHT_RARITY = 10             # aggregate rarity / protective alleles

# Complex-trait contribution thresholds
CONTRIBUTION_HIGH_THRESHOLD = 70
CONTRIBUTION_MODERATE_THRESHOLD = 40
CONTRIBUTION_LOW_THRESHOLD = 20

# Per-item scoring caps (complex mode)
COMPLEX_GWAS_HET_POINTS = 3            # points per lead SNP heterozygous genotype
COMPLEX_GWAS_HOM_POINTS = 6            # points per lead SNP homozygous genotype
COMPLEX_GWAS_MAX_POINTS = 25           # cap for GWAS dimension
COMPLEX_RARE_MISSENSE_BASE = 5         # base points for a rare missense in core gene
COMPLEX_RARE_LOF_BASE = 10             # base points for a rare LoF in core gene
COMPLEX_RARE_CLINVAR_LP_P_POINTS = 5   # bonus for ClinVar likely-pathogenic/pathogenic
COMPLEX_RARE_CLINVAR_BENIGN_PENALTY = -3  # penalty for ClinVar benign/likely-benign
COMPLEX_RARE_MAX_PER_VARIANT = 20      # per-variant cap in rare-functional dimension
COMPLEX_LIT_VARIANT_BONUS = 8          # variant-level literature support
COMPLEX_LIT_GENE_BONUS = 3             # gene-level literature support

# Disease mode defaults
DISEASE_MODE_MENDELIAN = "mendelian"
DISEASE_MODE_COMPLEX = "complex"
DISEASE_MODE_AUTO = "auto"
DEFAULT_DISEASE_MODE = DISEASE_MODE_AUTO

# APOE GRCh38 coordinates for Alzheimer risk assessment (1-based inclusive)
APOE_GRCH38 = {
    "chrom": "chr19",
    "start": 44908637,
    "end": 44912685,
    "rs7412": {"pos": 44908684, "ref": "T", "alt": "C"},  # ε2 risk / ε4 protective
    "rs429358": {"pos": 44908822, "ref": "C", "alt": "T"},  # ε4 risk
}

# Alzheimer disease gene tiers.
# Mendelian / high-effect familial genes get strict filtering and highest priority.
AD_MENDELIAN_CORE_GENES = frozenset({
    "APOE", "APP", "PSEN1", "PSEN2", "ABCA7", "TREM2", "SORL1", "ADAM10",
})
# GWAS / association loci are informative but should not relax AF thresholds,
# because many contain common HLA/polymorphic variants that explode noise.
AD_GWAS_LOCI_GENES = frozenset({
    # Established IGAP / AD GWAS loci
    "CLU", "PICALM", "BIN1", "CR1", "CD33", "MS4A4A", "EPHA1", "CD2AP",
    "INPP5D", "MEF2C", "HLA-DRB1", "PTK2B", "CASS4", "CELF1", "FERMT2",
    "SLC24A4", "RIN3", "DSG2", "PLCG2", "UNC5C", "AKAP9", "ADAMTS1",
    # Additional loci from PGS004863 (Sleiman et al. 2023, Kunkle 2019 GWAS)
    "ABCA1", "ABCA7", "ADAMTS1", "AKAP9", "ANK3", "ANKH", "APH1B", "APP",
    "BCKDK", "CTSB", "EML6", "FAM157C", "FAM171A2", "GRN", "HBEGF", "ICA1",
    "IDUA", "KLF7", "MAF", "MAPT", "MINDY2", "PLEKHA1", "PRKD3", "RASGEF1C",
    "RBCK1", "SHARPIN", "SNX1", "SPI1", "SPPL2A", "TPCN1", "TSPAN14", "WNT3",
})
# Backward-compatible union used for gene-set prioritization.
AD_CORE_GENES = AD_MENDELIAN_CORE_GENES | AD_GWAS_LOCI_GENES

# ClinVar phenotype matching keywords for common diseases
DISEASE_CLINVAR_KEYWORDS: dict[str, list[str]] = {
    "alzheimer disease": [
        "alzheimer", "alzheimer disease", "alzheimer's", "dementia",
        "cognitive decline", "memory loss", "presenile dementia",
    ],
    "parkinson disease": [
        "parkinson", "parkinson disease", "parkinson's", "parkinsonism",
    ],
    "cardiomyopathy": [
        "cardiomyopathy", "heart failure", "dilated cardiomyopathy",
        "hypertrophic cardiomyopathy", "arrhythmogenic cardiomyopathy",
    ],
    "hyperuricemia": [
        "hyperuricemia", "gout", "gouty arthritis", "urate",
        "uric acid", "nephrolithiasis", "urolithiasis", "kidney stone",
        "renal hypouricemia", "uric acid nephrolithiasis",
    ],
}

# Tier 3 denoising defaults
TIER3_MAX_GNOMAD_AF = 0.01          # drop common variants
TIER3_CORE_GENE_MAX_GNOMAD_AF = 0.05  # relaxed for AD core genes
TIER3_KEEP_IMPACTS = {"HIGH", "MODERATE"}
TIER3_MAX_VARIANTS_PER_GENE = 20    # collapse excess common/benign hits

# Disease reference cache: one directory per disease key.
# Stores aggregated core genes, known pathogenic variants, GWAS lead SNPs,
# and literature evidence. Generated on first use and reused unless refreshed.
DISEASE_REF_CACHE = WORKBUDDY_DATA / "drq_disease_refs"

# Built-in disease reference seeds (core genes + known GWAS/pathogenic loci).
# These are used to bootstrap the cache and prioritize the gene set.
DISEASE_BUILTIN_REFS: dict[str, dict] = {
    "alzheimer disease": {
        "mode": DISEASE_MODE_MENDELIAN,
        "core_genes": sorted(AD_MENDELIAN_CORE_GENES),
        "gwas_loci_genes": sorted(AD_GWAS_LOCI_GENES),
        "gwas_lead_snps": [
            # APOE and TREM2 familial/high-effect variants (not included in PGS004863)
            {"rsid": "rs429358", "chrom": "chr19", "pos": 44908822, "gene": "APOE", "note": "APOE ε4 tag"},
            {"rsid": "rs7412", "chrom": "chr19", "pos": 44908684, "gene": "APOE", "note": "APOE ε2 tag"},
            {"rsid": "rs75932628", "chrom": "chr6", "pos": 41161514, "gene": "TREM2", "note": "R47H risk variant"},
            # PGS004863 (Sleiman et al. 2023, Alzheimers Dement) AD PRS lead SNPs, GRCh38, n=74
            {"rsid": "rs11111149", "chrom": "chr12", "pos": 102072301, "gene": "ENSG00000257222", "note": "PGS004863 OR=1.27"},
            {"rsid": "rs6733839", "chrom": "chr2", "pos": 127135234, "gene": "BIN1", "note": "PGS004863 OR=1.19"},
            {"rsid": "rs745717", "chrom": "chr2", "pos": 127136321, "gene": "BIN1", "note": "PGS004863 OR=1.19"},
            {"rsid": "rs35114168", "chrom": "chr2", "pos": 127090354, "gene": "BIN1", "note": "PGS004863 OR=1.15"},
            {"rsid": "rs35220752", "chrom": "chr11", "pos": 85941686, "gene": "PICALM", "note": "PGS004863 OR=0.86"},
            {"rsid": "rs679515", "chrom": "chr1", "pos": 207577223, "gene": "CR1", "note": "PGS004863 OR=1.14"},
            {"rsid": "rs12998662", "chrom": "chr2", "pos": 207101667, "gene": "KLF7", "note": "PGS004863 OR=1.13"},
            {"rsid": "rs12618593", "chrom": "chr2", "pos": 127128795, "gene": "BIN1", "note": "PGS004863 OR=1.13"},
            {"rsid": "rs56407236", "chrom": "chr16", "pos": 90103687, "gene": "FAM157C", "note": "PGS004863 OR=1.13"},
            {"rsid": "rs117618017", "chrom": "chr15", "pos": 63277703, "gene": "APH1B", "note": "PGS004863 OR=1.13"},
            {"rsid": "rs1065712", "chrom": "chr8", "pos": 11844613, "gene": "CTSB", "note": "PGS004863 OR=1.13"},
            {"rsid": "rs9373079", "chrom": "chr6", "pos": 134068394, "gene": "ENSG00000287413", "note": "PGS004863 OR=1.12"},
            {"rsid": "rs12151021", "chrom": "chr19", "pos": 1050875, "gene": "ABCA7", "note": "PGS004863 OR=1.12"},
            {"rsid": "rs34173062", "chrom": "chr8", "pos": 144103704, "gene": "SHARPIN", "note": "PGS004863 OR=1.12"},
            {"rsid": "rs7564197", "chrom": "chr2", "pos": 127123792, "gene": "BIN1", "note": "PGS004863 OR=0.88"},
            {"rsid": "rs6014724", "chrom": "chr20", "pos": 56423488, "gene": "CASS4", "note": "PGS004863 OR=0.88"},
            {"rsid": "rs1060743", "chrom": "chr2", "pos": 127068957, "gene": "BIN1", "note": "PGS004863 OR=1.11"},
            {"rsid": "rs11525457", "chrom": "chr11", "pos": 60321826, "gene": "MS4A4A", "note": "PGS004863 OR=1.11"},
            {"rsid": "rs3851179", "chrom": "chr11", "pos": 86157598, "gene": "PICALM", "note": "PGS004863 OR=0.89"},
            {"rsid": "rs113706587", "chrom": "chr5", "pos": 180201150, "gene": "RASGEF1C", "note": "PGS004863 OR=1.1"},
            {"rsid": "rs11787077", "chrom": "chr8", "pos": 27607795, "gene": "CLU", "note": "PGS004863 OR=0.9"},
            {"rsid": "rs2972178", "chrom": "chr8", "pos": 1739731, "gene": "CLU", "note": "PGS004863 OR=0.9"},
            {"rsid": "rs17125924", "chrom": "chr14", "pos": 52924962, "gene": "FERMT2", "note": "PGS004863 OR=1.09"},
            {"rsid": "rs112403360", "chrom": "chr5", "pos": 14724304, "gene": "ANKH", "note": "PGS004863 OR=1.09"},
            {"rsid": "rs1582763", "chrom": "chr11", "pos": 60254475, "gene": "MS4A4A", "note": "PGS004863 OR=0.91"},
            {"rsid": "rs2632516", "chrom": "chr17", "pos": 58331728, "gene": "MAPT", "note": "PGS004863 OR=0.91"},
            {"rsid": "rs543928", "chrom": "chr6", "pos": 32317471, "gene": "HLA-DRB1", "note": "PGS004863 OR=1.08"},
            {"rsid": "rs598561", "chrom": "chr11", "pos": 85941783, "gene": "PICALM", "note": "PGS004863 OR=1.08"},
            {"rsid": "rs11168036", "chrom": "chr5", "pos": 140327854, "gene": "HBEGF", "note": "PGS004863 OR=1.08"},
            {"rsid": "rs7225151", "chrom": "chr17", "pos": 5233752, "gene": "ZNF594-DT", "note": "PGS004863 OR=1.08"},
            {"rsid": "rs6489896", "chrom": "chr12", "pos": 113281983, "gene": "TPCN1", "note": "PGS004863 OR=1.08"},
            {"rsid": "rs62374257", "chrom": "chr5", "pos": 86927378, "gene": "ENSG00000249061", "note": "PGS004863 OR=1.08"},
            {"rsid": "rs12358692", "chrom": "chr10", "pos": 11679103, "gene": "ENSG00000271046", "note": "PGS004863 OR=0.92"},
            {"rsid": "rs6605556", "chrom": "chr6", "pos": 32615322, "gene": "HLA-DRB1", "note": "PGS004863 OR=0.92"},
            {"rsid": "rs6586028", "chrom": "chr10", "pos": 80494228, "gene": "TSPAN14", "note": "PGS004863 OR=0.92"},
            {"rsid": "rs10933431", "chrom": "chr2", "pos": 233117202, "gene": "INPP5D", "note": "PGS004863 OR=0.92"},
            {"rsid": "rs73223431", "chrom": "chr8", "pos": 27362470, "gene": "PTK2B", "note": "PGS004863 OR=1.07"},
            {"rsid": "rs13391802", "chrom": "chr2", "pos": 54918863, "gene": "EML6", "note": "PGS004863 OR=1.07"},
            {"rsid": "rs7157106", "chrom": "chr14", "pos": 105761758, "gene": "IGHG3", "note": "PGS004863 OR=1.07"},
            {"rsid": "rs587709", "chrom": "chr19", "pos": 54267597, "gene": "SPI1", "note": "PGS004863 OR=1.07"},
            {"rsid": "rs4985556", "chrom": "chr16", "pos": 70660097, "gene": "PLCG2", "note": "PGS004863 OR=1.07"},
            {"rsid": "rs1800978", "chrom": "chr9", "pos": 104903697, "gene": "ABCA1", "note": "PGS004863 OR=1.07"},
            {"rsid": "rs10952097", "chrom": "chr7", "pos": 8204382, "gene": "ICA1", "note": "PGS004863 OR=1.07"},
            {"rsid": "rs602602", "chrom": "chr15", "pos": 58764824, "gene": "MINDY2", "note": "PGS004863 OR=0.93"},
            {"rsid": "rs12590654", "chrom": "chr14", "pos": 92472511, "gene": "SLC24A4", "note": "PGS004863 OR=0.93"},
            {"rsid": "rs11771145", "chrom": "chr7", "pos": 143413669, "gene": "EPHA1-AS1", "note": "PGS004863 OR=0.93"},
            {"rsid": "rs10131280", "chrom": "chr14", "pos": 106665591, "gene": "IGHV3-65", "note": "PGS004863 OR=0.93"},
            {"rsid": "rs7912495", "chrom": "chr10", "pos": 11676714, "gene": "FAM171A2", "note": "PGS004863 OR=1.06"},
            {"rsid": "rs7767350", "chrom": "chr6", "pos": 47517390, "gene": "CD2AP", "note": "PGS004863 OR=1.06"},
            {"rsid": "rs76928645", "chrom": "chr7", "pos": 54873635, "gene": "AKAP9", "note": "PGS004863 OR=0.94"},
            {"rsid": "rs7384878", "chrom": "chr7", "pos": 100334426, "gene": "EPHA1", "note": "PGS004863 OR=0.94"},
            {"rsid": "rs72777026", "chrom": "chr2", "pos": 9558882, "gene": "ENSG00000271855", "note": "PGS004863 OR=1.06"},
            {"rsid": "rs6943429", "chrom": "chr7", "pos": 7817263, "gene": "CD2AP", "note": "PGS004863 OR=1.06"},
            {"rsid": "rs5848", "chrom": "chr17", "pos": 44352876, "gene": "GRN", "note": "PGS004863 OR=1.06"},
            {"rsid": "rs4277405", "chrom": "chr17", "pos": 63471557, "gene": "PPIAP55", "note": "PGS004863 OR=0.94"},
            {"rsid": "rs2245466", "chrom": "chr4", "pos": 40197226, "gene": "UNC5C", "note": "PGS004863 OR=1.06"},
            {"rsid": "rs199515", "chrom": "chr17", "pos": 46779275, "gene": "WNT3", "note": "PGS004863 OR=0.94"},
            {"rsid": "rs17020490", "chrom": "chr2", "pos": 37304796, "gene": "PRKD3", "note": "PGS004863 OR=1.06"},
            {"rsid": "rs1358782", "chrom": "chr20", "pos": 413334, "gene": "RBCK1", "note": "PGS004863 OR=0.94"},
            {"rsid": "rs12446759", "chrom": "chr16", "pos": 81739398, "gene": "FAM157C", "note": "PGS004863 OR=0.94"},
            {"rsid": "rs10437655", "chrom": "chr11", "pos": 47370397, "gene": "SPI1", "note": "PGS004863 OR=1.06"},
            {"rsid": "rs9304690", "chrom": "chr19", "pos": 49950060, "gene": "ENSG00000269179", "note": "PGS004863 OR=1.05"},
            {"rsid": "rs6966331", "chrom": "chr7", "pos": 37844191, "gene": "ENSG00000290149", "note": "PGS004863 OR=0.95"},
            {"rsid": "rs6846529", "chrom": "chr4", "pos": 11023507, "gene": "LINC02498", "note": "PGS004863 OR=1.05"},
            {"rsid": "rs450674", "chrom": "chr16", "pos": 79574511, "gene": "MAF", "note": "PGS004863 OR=0.95"},
            {"rsid": "rs3848143", "chrom": "chr15", "pos": 64131307, "gene": "SNX1", "note": "PGS004863 OR=1.05"},
            {"rsid": "rs3822030", "chrom": "chr4", "pos": 993555, "gene": "IDUA", "note": "PGS004863 OR=0.95"},
            {"rsid": "rs2830489", "chrom": "chr21", "pos": 26775872, "gene": "ADAMTS1", "note": "PGS004863 OR=0.95"},
            {"rsid": "rs2154481", "chrom": "chr21", "pos": 26101558, "gene": "APP", "note": "PGS004863 OR=0.95"},
            {"rsid": "rs889555", "chrom": "chr16", "pos": 31111250, "gene": "BCKDK", "note": "PGS004863 OR=0.96"},
            {"rsid": "rs8025980", "chrom": "chr15", "pos": 50701814, "gene": "SPPL2A", "note": "PGS004863 OR=0.96"},
            {"rsid": "rs7908662", "chrom": "chr10", "pos": 122413396, "gene": "PLEKHA1", "note": "PGS004863 OR=0.96"},
            {"rsid": "rs7401792", "chrom": "chr14", "pos": 92464917, "gene": "SLC24A4", "note": "PGS004863 OR=1.04"},
            {"rsid": "rs7068231", "chrom": "chr10", "pos": 60025170, "gene": "ANK3", "note": "PGS004863 OR=0.96"},
        ],
        "key_literature": [
            {"pmid": "15703407", "title": "A high-density whole-genome association study reveals ...", "genes": ["CLU", "PICALM"], "note": "Lambert et al., Nat Genet 2009"},
            {"pmid": "21529783", "title": "Common variants at MS4A4/MS4A6E, CD2AP, CD33 and EPHA1 are associated with Alzheimer's disease.", "genes": ["MS4A4A", "CD2AP", "CD33", "EPHA1"], "note": "Hooli et al. / Naj et al., Nat Genet 2011"},
            {"pmid": "23184178", "title": "TREM2 variants in Alzheimer's disease.", "genes": ["TREM2"], "note": "Guerreiro et al., N Engl J Med 2013"},
            {"pmid": "12499452", "title": "SORL1 is genetically associated with late-onset Alzheimer's disease.", "genes": ["SORL1"], "note": "Rogaeva et al., Nat Genet 2007"},
            {"pmid": "19914914", "title": "Meta-analysis confirms CR1, CLU, and PICALM as Alzheimer's disease risk loci.", "genes": ["CR1", "CLU", "PICALM"], "note": "Jun et al., Alzheimers Dement 2010"},
            {"pmid": "24162737", "title": "GWAS of 74,046 individuals identifies 11 new susceptibility loci for Alzheimer's disease.", "genes": ["HLA-DRB1", "PTK2B", "CASS4", "CELF1", "FERMT2", "SLC24A4", "RIN3", "DSG2", "PLCG2", "UNC5C"], "note": "Lambert et al., Nat Genet 2013 (IGAP stage 1)"},
            {"pmid": "30617256", "title": "New insights into the genetic etiology of Alzheimer's disease.", "genes": ["ABCA7", "INPP5D", "MEF2C", "AKAP9", "ADAMTS1"], "note": "Jansen et al., Nat Genet 2019 (stage 2)"},
            {"pmid": "8290044", "title": "Gene dose of apolipoprotein E type 4 allele and the risk of Alzheimer's disease.", "genes": ["APOE"], "note": "Corder et al., Science 1993"},
            {"pmid": "25700176", "title": "Rare coding variants in PLCG2, ABI3, and TREM2 implicate microglial-mediated innate immunity.", "genes": ["PLCG2"], "note": "Sims et al., Nat Genet 2017"},
            {"pmid": "37450379", "title": "Trans-ethnic genomic informed risk assessment for Alzheimer's disease: An International Hundred K+ Cohorts Consortium study.", "genes": ["BIN1", "PICALM", "CLU", "CR1", "ABCA7", "CD2AP", "INPP5D", "MEF2C", "HLA-DRB1", "PTK2B", "CASS4", "FERMT2", "SLC24A4", "RIN3", "PLCG2", "UNC5C", "AKAP9", "ADAMTS1", "TSPAN14", "GRN", "SPI1"], "note": "Sleiman et al., Alzheimers Dement 2023 (PGS004863)"},
        ],
    },
    "parkinson disease": {
        "mode": DISEASE_MODE_MENDELIAN,
        "core_genes": ["SNCA", "LRRK2", "PARK7", "PINK1", "PRKN", "GBA1", "VPS35", "DJ1", "UCHL1", "ATP13A2", "PLA2G6", "FBXO7", "DNAJC6", "SYNJ1", "VPS13C", "CHCHD2", "GBA", "MAPT"],
        "gwas_lead_snps": [
            {"rsid": "rs34637584", "chrom": "chr12", "pos": 40340400, "gene": "LRRK2", "note": "G2019S risk"},
            {"rsid": "rs11931074", "chrom": "chr4", "pos": 90645145, "gene": "SNCA", "note": "SNCA locus"},
            {"rsid": "rs17649553", "chrom": "chr1", "pos": 155209060, "gene": "GBA1", "note": "GBA1 locus"},
            {"rsid": "rs393152", "chrom": "chr17", "pos": 45912131, "gene": "MAPT", "note": "MAPT locus"},
        ],
        "key_literature": [
            {"pmid": "16406073", "title": "PINK1 mutations are associated with sporadic early-onset parkinsonism.", "genes": ["PINK1"], "note": "Valente et al."},
            {"pmid": "18509534", "title": "Mutations in LRRK2 cause autosomal-dominant parkinsonism ...", "genes": ["LRRK2"], "note": "Zimprich et al."},
        ],
    },
    "adult vision disorders": {
        "mode": DISEASE_MODE_COMPLEX,
        "core_genes": [
            # Age-related macular degeneration (AMD)
            "CFH", "ARMS2", "HTRA1", "C2", "CFB", "C3", "CFI", "APOE", "CETP", "LIPC",
            "TIMP3", "VEGFA", "COL8A1", "RAD51B", "TNFRSF10A", "ARL6IP6", "PILRB", "PILRA",
            # 2024 multi-ancestry AMD GWAS novel loci (PMC12817227)
            "CD46", "CD55", "PLTP", "MMP9", "CFD", "RRAS", "SERPINA1", "TGFB1", "SMAD3",
            "ADAM19", "TYR", "HERC2", "OCA2", "TRPM1", "RASIP1", "IGFBP7", "PDGFB", "MYO1E",
            "EXOC3L2", "CSK", "ULK3", "ZBTB38", "ZNF385B", "CAND2", "TMEM40", "CHD9", "LBP",
            "HSDL2", "AFF1", "ACAA2", "LIPG", "ME3", "RLBP1", "CLUL1", "LRP2", "RP1L1", "C5",
            "HLA-DQB1", "HLA-DRB1",
            # Primary open-angle glaucoma (POAG) / normal-tension glaucoma (NTG)
            "MYOC", "OPTN", "WDR36", "TBK1", "CYP1B1", "FOXC1", "PITX2", "PAX6", "LTBP2",
            "TEK", "ANGPT1", "THBS1", "EFEMP1", "LOXL1", "CDKN2B-AS1", "TMCO1", "ABCA1",
            "CAV1", "CAV2", "SIX1", "SIX6", "GAS7", "FOXO1", "ATOH7",
            # 2024/2025 POAG novel loci (All of Us / multi-ancestry)
            "TUT4", "RYK", "MOXD1", "UBAP2", "PRRX1", "TSPAN17", "SLC16A7", "LINC02388",
            "FAM135B", "LINC00871", "GATA5", "SGCZ", "MAFTRR", "GGT7", "TRPM2",
            # High myopia / refractive error
            "PAX6", "MYP1", "MYP2", "MYP3", "MYP4", "MYP5", "MYP6", "MYP9", "MYP10",
            "MYP11", "MYP12", "MYP15", "MYP16", "ZNF644", "SCO2", "CCND1", "BICC1",
            "CD55", "KCNQ5", "GJD2", "RBFOX1", "RDH5", "BMP3", "PRSS56", "MPPED2", "CYP26A1",
            # Common retinal/photoreceptor genes with adult-onset phenotypes
            "ABCA4", "EYS", "USH2A", "RPGR", "RHO", "PRPH2", "CHM",
        ],
        "gwas_loci_genes": [
            "CFH", "ARMS2", "HTRA1", "C3", "CFB", "C2", "APOE", "VEGFA",
            "CD46", "CD55", "PLTP", "CFD", "RRAS", "SERPINA1", "TGFB1", "SMAD3", "ADAM19",
            "TYR", "HERC2", "OCA2", "TRPM1", "RASIP1", "IGFBP7", "PDGFB", "MYO1E",
            "CSK", "ULK3", "ZBTB38", "ZNF385B", "CHD9", "LBP", "HSDL2", "ACAA2", "LIPG",
            "RLBP1", "CLUL1", "LRP2", "RP1L1", "C5",
            "MYOC", "CDKN2B-AS1", "TMCO1", "SIX1", "SIX6", "ABCA1", "CAV1", "CAV2",
            "TUT4", "RYK", "MOXD1", "UBAP2", "PRRX1", "TSPAN17", "SLC16A7", "FAM135B",
            "GATA5", "SGCZ", "MAFTRR", "GGT7", "TRPM2",
            "PAX6", "GJD2", "RBFOX1", "KCNQ5", "BMP3", "MPPED2", "CYP26A1",
        ],
        "gwas_lead_snps": [
            # AMD
            {"rsid": "rs1061170", "chrom": "chr1", "pos": 196659237, "gene": "CFH", "note": "CFH Y402H, major AMD risk variant"},
            {"rsid": "rs10490924", "chrom": "chr10", "pos": 124204438, "gene": "ARMS2", "note": "ARMS2 A69S / HTRA1 locus"},
            {"rsid": "rs2230199", "chrom": "chr19", "pos": 6679385, "gene": "C3", "note": "C3 R102G"},
            {"rsid": "rs10922109", "chrom": "chr13", "pos": 110430300, "gene": "CFI", "note": "CFI locus"},
            {"rsid": "rs429358", "chrom": "chr19", "pos": 44908822, "gene": "APOE", "note": "APOE ε4, AMD risk"},
            {"rsid": "rs7412", "chrom": "chr19", "pos": 44908684, "gene": "APOE", "note": "APOE ε2"},
            {"rsid": "rs13278062", "chrom": "chr8", "pos": 23225418, "gene": "TNFRSF10A", "note": "TNFRSF10A AMD locus"},
            {"rsid": "rs5754227", "chrom": "chr22", "pos": 33157620, "gene": "SYN3", "note": "SYN3/TIMP3 AMD locus"},
            # 2024 multi-ancestry AMD GWAS (PMC12817227)
            {"rsid": "rs2724360", "chrom": "chr1", "pos": 207769813, "gene": "CD46/CD55", "note": "Novel complement-regulation locus"},
            {"rsid": "rs17447545", "chrom": "chr20", "pos": 45918429, "gene": "PLTP/MMP9", "note": "Novel lipid/angiogenesis locus"},
            {"rsid": "rs35186399", "chrom": "chr19", "pos": 860766, "gene": "CFD", "note": "Rare protective missense p.Glu69Lys"},
            {"rsid": "rs61760904", "chrom": "chr19", "pos": 49636675, "gene": "RRAS", "note": "Rare protective missense p.Asp133Asn"},
            {"rsid": "rs28929474", "chrom": "chr14", "pos": 94378610, "gene": "SERPINA1", "note": "PiZ, protective"},
            {"rsid": "rs17580", "chrom": "chr14", "pos": 94380925, "gene": "SERPINA1", "note": "PiS, risk"},
            {"rsid": "rs3844313", "chrom": "chr6", "pos": 32667852, "gene": "HLA-DQB1", "note": "African-ancestry MHC class II signal"},
            {"rsid": "rs28383172", "chrom": "chr6", "pos": 32598202, "gene": "HLA-DRB1", "note": "Risk haplotype tag"},
            {"rsid": "rs7775228", "chrom": "chr6", "pos": 32690302, "gene": "HLA-DQB1", "note": "MHC class II risk haplotype tag"},
            # POAG / glaucoma
            {"rsid": "rs2157719", "chrom": "chr9", "pos": 22003367, "gene": "CDKN2B-AS1", "note": "POAG GWAS locus"},
            {"rsid": "rs7555523", "chrom": "chr14", "pos": 60788658, "gene": "SIX1", "note": "SIX1/SIX6 optic nerve locus"},
            {"rsid": "rs1900004", "chrom": "chr9", "pos": 107560750, "gene": "ABCA1", "note": "ABCA1 POAG locus"},
            {"rsid": "rs7081455", "chrom": "chr1", "pos": 71500609, "gene": "PLXDC2", "note": "PLXDC2/TMCO1 POAG locus"},
            {"rsid": "rs1052990", "chrom": "chr7", "pos": 115993814, "gene": "CAV1", "note": "CAV1/CAV2 POAG locus"},
            {"rsid": "rs11656696", "chrom": "chr1", "pos": 171625297, "gene": "GAS7", "note": "GAS7 POAG locus"},
            {"rsid": "rs560766", "chrom": "chr1", "pos": 47299151, "gene": "ATOH7", "note": "ATOH7 optic disc locus"},
            # High myopia / refractive error
            {"rsid": "rs4373767", "chrom": "chr1", "pos": 84464929, "gene": "GJD2", "note": "GJD2 refractive error locus"},
            {"rsid": "rs524952", "chrom": "chr11", "pos": 30309797, "gene": "MPPED2", "note": "MPPED2 myopia locus"},
            {"rsid": "rs6888455", "chrom": "chr5", "pos": 77268628, "gene": "KCNQ5", "note": "KCNQ5 refractive error locus"},
            {"rsid": "rs634990", "chrom": "chr15", "pos": 67410633, "gene": "CYP26A1", "note": "CYP26A1 myopia locus"},
            {"rsid": "rs4134557", "chrom": "chr2", "pos": 31995, "gene": "RBFOX1", "note": "RBFOX1 refractive error locus"},
            {"rsid": "rs28500726", "chrom": "chr12", "pos": 44334700, "gene": "BMP3", "note": "BMP3 myopia locus"},
        ],
        "key_literature": [
            {"pmid": "16174643", "title": "Complement factor H polymorphism in age-related macular degeneration.", "genes": ["CFH"], "note": "Edwards et al., Science 2005"},
            {"pmid": "17486057", "title": "Common variation in three genes, including a noncoding variant in CFH, strongly influences risk of age-related macular degeneration.", "genes": ["CFH", "ARMS2", "HTRA1"], "note": "Maller et al., Nat Genet 2007"},
            {"pmid": "20385819", "title": "Common variants near FRK/COL10A1 and VEGFA are associated with advanced age-related macular degeneration.", "genes": ["VEGFA", "COL8A1"], "note": "Chen et al., Hum Mol Genet 2010"},
            {"pmid": "26691988", "title": "Seven new loci associated with age-related macular degeneration.", "genes": ["CFI", "CETP", "LIPC", "TIMP3", "RAD51B", "TNFRSF10A"], "note": "Fritsche et al., Nat Genet 2016"},
            {"pmid": "39590698", "title": "Multi-ancestry genome-wide association study of age-related macular degeneration identifies 63 loci, including 30 novel.", "genes": ["CD46", "CD55", "PLTP", "CFD", "RRAS", "SERPINA1", "TGFB1", "SMAD3", "ADAM19", "TYR", "HERC2", "OCA2", "TRPM1", "RASIP1", "IGFBP7", "PDGFB", "MYO1E", "EXOC3L2", "CSK", "ULK3", "ZBTB38", "ZNF385B", "CAND2", "TMEM40", "CHD9", "LBP", "HSDL2", "AFF1", "ACAA2", "LIPG", "ME3", "RLBP1", "CLUL1", "LRP2", "RP1L1", "C5", "HLA-DQB1", "HLA-DRB1"], "note": "PMID:39590698 / PMC12817227, 2024"},
            {"pmid": "21685912", "title": "Common variants at 9p21 and 8q22 are associated with increased susceptibility to optic nerve degeneration in glaucoma.", "genes": ["CDKN2B-AS1"], "note": "Wiggs et al., PLoS Genet 2012"},
            {"pmid": "26691985", "title": "Genome-wide association analyses identify 139 new loci for primary open-angle glaucoma.", "genes": ["MYOC", "ABCA1", "SIX1", "SIX6", "GAS7", "CAV1", "CAV2", "TMCO1"], "note": "Gharahkhani et al., Nat Genet 2024 / IGAG"},
            {"pmid": "40049121", "title": "Multi-ancestry genome-wide association study in All of Us for primary open-angle glaucoma.", "genes": ["TUT4", "RYK", "MOXD1", "UBAP2", "TSPAN17", "SLC16A7", "FAM135B", "GATA5", "SGCZ", "MAFTRR", "GGT7", "TRPM2"], "note": "PMID:40049121 / PMC13129092, 2025"},
            {"pmid": "20835237", "title": "Common variants on chromosome 5p12 and 1q41 confer susceptibility to high myopia.", "genes": ["MYP10", "MYP15"], "note": "Shi et al., Nat Commun 2011"},
            {"pmid": "20835238", "title": "A genome-wide association study identifies a susceptibility locus for refractive errors and myopia at 15q14.", "genes": ["GJD2"], "note": "Solouki et al., Nat Genet 2010"},
        ],
    },
    "hyperuricemia": {
        "mode": DISEASE_MODE_COMPLEX,
        "core_genes": [
            # Urate transporters and renal handling
            "SLC2A9", "SLC22A12", "SLC22A11", "ABCG2", "SLC17A1", "SLC17A3", "PDZK1",
            # Purine metabolism
            "HPRT1", "PRPS1", "XDH", "MOCOS", "AOX1",
            # Renal tubular / electrolyte genes that alter urate handling
            "UMOD", "REN", "HNF1B", "CLDN16", "CLDN19",
            # Metabolic syndrome / lipids / ketone bodies influencing urate
            "GCKR", "FADS1", "PPARG", "APOB", "APOC3", "LPL",
            "HMGCL", "HMGCS2", "OXCT1", "BDH1", "ACAT1", "ACAT2", "HMGCS1",
            # Additional GWAS loci
            "ALDH16A1", "RFX3", "A1CF", "INHBC", "SF1", "TRIM46",
        ],
        "gwas_loci_genes": [
            "SLC2A9", "SLC22A12", "ABCG2", "SLC17A1", "SLC17A3", "PDZK1",
            "GCKR", "RFX3", "A1CF", "INHBC", "SF1", "TRIM46", "ALDH16A1",
        ],
        "gwas_lead_snps": [
            {"rsid": "rs12498742", "chrom": "chr4", "pos": 9942428, "gene": "SLC2A9", "note": "SLC2A9/GLUT9 major urate transporter locus"},
            {"rsid": "rs1014290", "chrom": "chr4", "pos": 10000237, "gene": "SLC2A9", "note": "SLC2A9 secondary signal"},
            {"rsid": "rs6855911", "chrom": "chr4", "pos": 9934286, "gene": "SLC2A9", "note": "SLC2A9 third signal"},
            {"rsid": "rs505802", "chrom": "chr11", "pos": 64589600, "gene": "SLC22A12", "note": "SLC22A12/URAT1 locus"},
            {"rsid": "rs7932775", "chrom": "chr11", "pos": 64600390, "gene": "SLC22A12", "note": "URAT1 secondary signal"},
            {"rsid": "rs12129861", "chrom": "chr1", "pos": 145709377, "gene": "PDZK1", "note": "PDZK1 scaffolding protein, renal urate handling"},
            {"rsid": "rs2231142", "chrom": "chr4", "pos": 88131171, "gene": "ABCG2", "note": "ABCG2 Q141K, gut/renal urate excretion"},
            {"rsid": "rs3114018", "chrom": "chr4", "pos": 88143429, "gene": "ABCG2", "note": "ABCG2 secondary signal"},
            {"rsid": "rs1165196", "chrom": "chr6", "pos": 25812922, "gene": "SLC17A1", "note": "SLC17A1/NPT1 renal urate secretion"},
            {"rsid": "rs1183201", "chrom": "chr6", "pos": 25823216, "gene": "SLC17A1", "note": "SLC17A1 secondary signal"},
            {"rsid": "rs1165205", "chrom": "chr6", "pos": 25870314, "gene": "SLC17A3", "note": "SLC17A3/NPT4 renal urate secretion"},
            {"rsid": "rs780094", "chrom": "chr2", "pos": 27518370, "gene": "GCKR", "note": "GCKR glucokinase regulator, glucose/urate metabolism"},
            {"rsid": "rs1260326", "chrom": "chr2", "pos": 27508073, "gene": "GCKR", "note": "GCKR P446L, triglyceride/urate link"},
            {"rsid": "rs10821905", "chrom": "chr10", "pos": 50886333, "gene": "RFX3", "note": "RFX3 pancreatic islet / urate GWAS locus"},
            {"rsid": "rs2078267", "chrom": "chr11", "pos": 64566642, "gene": "A1CF", "note": "A1CF / SLC22A12 region urate locus"},
            {"rsid": "rs12356193", "chrom": "chr10", "pos": 59653595, "gene": "INHBC", "note": "INHBC urate GWAS locus"},
            {"rsid": "rs12801356", "chrom": "chr11", "pos": 122830307, "gene": "SF1", "note": "SF1 urate GWAS locus"},
            {"rsid": "rs16998073", "chrom": "chr4", "pos": 80263187, "gene": "ALDH16A1", "note": "ALDH16A1 urate GWAS locus"},
        ],
        "key_literature": [
            {"pmid": "18179892", "title": "SLC2A9 is a newly identified urate transporter influencing serum urate concentration, urate excretion and gout.", "genes": ["SLC2A9"], "note": "Vitart et al., Nat Genet 2008"},
            {"pmid": "19503597", "title": "Common dysfunctional variants in ABCG2 are a major cause of early-onset gout.", "genes": ["ABCG2"], "note": "Matsuo et al., Sci Rep 2011"},
            {"pmid": "23263486", "title": "Genome-wide association analyses identify 18 new loci associated with serum urate concentrations.", "genes": ["SLC2A9", "SLC22A12", "ABCG2", "SLC17A1", "SLC17A3", "GCKR", "RFX3", "PDZK1"], "note": "Köttgen et al., Nat Genet 2013"},
            {"pmid": "23785249", "title": "Genome-wide association study of clinically defined gout identifies multiple risk loci and its association with serum urate levels.", "genes": ["SLC2A9", "ABCG2", "SLC22A12", "PDZK1", "GCKR"], "note": "Okada et al., Ann Rheum Dis 2012"},
            {"pmid": "19506252", "title": "Genetic variants in SLC2A9, SLC22A12 and ABCG2 associate with diuretic-induced gout.", "genes": ["SLC2A9", "SLC22A12", "ABCG2"], "note": "McKeigue et al., Pharmacogenomics 2009"},
            {"pmid": "26467716", "title": "An update on the genetics of hyperuricemia and gout.", "genes": ["SLC2A9", "SLC22A12", "ABCG2", "SLC17A1", "SLC17A3", "GCKR"], "note": "Merriman & Dalbeth, Curr Opin Rheumatol 2015"},
            {"pmid": "29456098", "title": "A genome-wide association study of gout in multiple ethnic groups.", "genes": ["SLC2A9", "ABCG2", "SLC22A12", "PDZK1", "GCKR", "ALDH16A1"], "note": "Tin et al., Arthritis Rheumatol 2018"},
        ],
    },
}

# Diseases for which a built-in reference is available.
BUILTIN_DISEASE_KEYS = set(DISEASE_BUILTIN_REFS.keys())

# Keywords used in auto mode to classify a disease as complex-trait.
COMPLEX_DISEASE_KEYWORDS = [
    "hyperuricemia", "gout", "urate",
    "diabetes", "type 2 diabetes", "t2d", "glucose", "insulin resistance",
    "obesity", "bmi", "body mass index", "adiposity",
    "dyslipidemia", "hyperlipidemia", "cholesterol", "triglyceride", "ldl", "hdl",
    "hypertension", "blood pressure",
    "metabolic syndrome", "metabolic",
    "non-alcoholic fatty liver", "nafld", "nash",
    "osteoporosis", "bone mineral density", "bmd",
    "height", "stature", "growth",
    "longevity", "lifespan",
    "intelligence", "iq", "cognitive ability",
    "personality", "depression", "anxiety",
    "autoimmune", "rheumatoid arthritis", "lupus", "sle", "crohn", "ulcerative colitis",
    "asthma", "allergy", "atopic",
    "vision", "macular", "age-related macular degeneration", "amd",
    "glaucoma", "optic nerve", "intraocular pressure",
    "myopia", "refractive error", "high myopia", "cataract", "retinopathy",
]


def resolve_disease_mode(disease_name: str, requested_mode: str = DISEASE_MODE_AUTO) -> str:
    """Resolve the disease mode for scoring/reporting.

    If the user explicitly requests mendelian or complex, return that.
    In auto mode, use built-in template mode when available; otherwise fall back
    to keyword heuristics for common complex traits.
    """
    if requested_mode in (DISEASE_MODE_MENDELIAN, DISEASE_MODE_COMPLEX):
        return requested_mode
    norm = disease_name.lower().strip()
    builtin = DISEASE_BUILTIN_REFS.get(norm)
    if builtin:
        return builtin.get("mode", DISEASE_MODE_MENDELIAN)
    # Try substring match against built-in keys
    for key, ref in DISEASE_BUILTIN_REFS.items():
        if key in norm or norm in key:
            return ref.get("mode", DISEASE_MODE_MENDELIAN)
    # Keyword heuristic
    for kw in COMPLEX_DISEASE_KEYWORDS:
        if kw in norm:
            return DISEASE_MODE_COMPLEX
    return DISEASE_MODE_MENDELIAN
