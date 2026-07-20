"""Body mass index (BMI) continuous-trait prediction template.

Source: Locke et al. 2015 (Nature) "Genetic studies of body mass index yield new insights for obesity biology"
        GIANT consortium BMI EUR ancestry summary statistics (SNP_gwas_mc_merge_nogc.tbl.uniq.gz).
        Wen et al. 2012 (Nat Genet) for East Asian BMI loci (CDKAL1, PCSK1, SEC16B, etc.).
        Banerjee & Girirajan 2025 (Nat Commun) for cross-ancestry rare-variant gene discovery.

Notes:
- This template predicts BMI as a continuous polygenic trait, not a binary obesity case/control status.
- Mendelian obesity genes (MC4R, LEP, LEPR, POMC, etc.) are included as high-impact contributors.
- Core GWAS loci (FTO, MC4R, BDNF, TMEM18, GNPDA2, NEGR1, MTCH2) are placed in prs_variants_high.
- Remaining validated BMI GWAS lead SNPs are placed in gwas_lead_snps.
- East Asian-specific loci (CDKAL1 rs2206734, PCSK1 rs261967, SEC16B rs2605100) are annotated
  with approximate EAS effect directions; EUR/EAS beta differences are noted per variant.
- Cross-ancestry genes from Banerjee & Girirajan 2025 (YLPM1, GIGYF1, SLC5A3, GRM7, BSN)
  are included in gene_set with low contribution scores.
- This template is independent from the T2D template to avoid double-counting FTO/T2D overlap.
"""

TEMPLATE_NAME = "body mass index"

TEMPLATE = {
    "aliases": [
        "body mass index",
        "bmi",
        "体重指数",
        "obesity",
        "肥胖",
        "肥胖倾向"
    ],
    "mode": "complex",
    "gene_set": [
        {
            "gene": "MC4R",
            "tier": "mendelian_mod",
            "contribution_score": 0.9,
            "penetrance": "high",
            "penetrance_score": 0.65,
            "evidence": "monogenic",
            "note": "黑素皮质素4受体;下丘脑食欲/能量稳态核心;变异导致显性肥胖综合征",
            "phenotype_assoc": "MC4R相关肥胖综合征;常显;外显率可变(~30-80%);严重儿童期肥胖",
            "key_domains": "7次跨膜GPCR结构域(aa1-333);N端胞外域(aa1-40);配体结合域(aa42-130);G蛋白偶联域(aa220-333)",
            "clingen_validity": "Definitive"
        },
        {
            "gene": "LEP",
            "tier": "mendelian_high",
            "contribution_score": 1.0,
            "penetrance": ">0.95",
            "penetrance_score": 0.95,
            "evidence": "monogenic",
            "note": "瘦素;脂肪细胞分泌;JAK-STAT信号;LOF导致严重早发肥胖",
            "phenotype_assoc": "先天性瘦素缺乏;常隐;严重早发肥胖、低促性腺激素性性腺功能减退",
            "key_domains": "信号肽(aa1-21);核心瘦素(aa22-167);二硫键(aa96-cys117)",
            "clingen_validity": "Definitive"
        },
        {
            "gene": "LEPR",
            "tier": "mendelian_high",
            "contribution_score": 1.0,
            "penetrance": ">0.95",
            "penetrance_score": 0.95,
            "evidence": "monogenic",
            "note": "瘦素受体;下丘脑表达;JAK-STAT/ERK信号;LOF导致瘦素抵抗型肥胖",
            "phenotype_assoc": "先天性瘦素受体缺乏;常隐;严重早发肥胖、促甲状腺激素/生长激素轴异常",
            "key_domains": "胞外瘦素结合域(aa428-635);纤维连接蛋白III型重复域;JAK2结合BOX motifs",
            "clingen_validity": "Definitive"
        },
        {
            "gene": "POMC",
            "tier": "mendelian_high",
            "contribution_score": 1.0,
            "penetrance": ">0.95",
            "penetrance_score": 0.95,
            "evidence": "monogenic",
            "note": "阿黑皮素原;垂体/下丘脑;裂解为α-MSH/ACTH/β-内啡肽;LOF导致食欲亢进",
            "phenotype_assoc": "POMC缺陷肥胖;常隐;肾上腺功能不全、红发、早发肥胖",
            "key_domains": "信号肽(aa1-26);ACTH(aa138-176);α-MSH(aa138-150);β-MSH(aa205-235);β-内啡肽",
            "clingen_validity": "Definitive"
        },
        {
            "gene": "SIM1",
            "tier": "mendelian_mod",
            "contribution_score": 0.7,
            "penetrance": "moderate",
            "penetrance_score": 0.4,
            "evidence": "monogenic",
            "note": "SIM1转录因子;下丘脑室旁核发育;调控MC4R表达;缺失导致肥胖综合征",
            "phenotype_assoc": "SIM1相关肥胖综合征/Prader-Willi样综合征;常显;食欲亢进、发育迟缓",
            "key_domains": "bHLH域(aa85-140);PAS域(aa270-360);转录激活域;核定位信号",
            "clingen_validity": "Strong"
        },
        {
            "gene": "BDNF",
            "tier": "mendelian_mod",
            "contribution_score": 0.7,
            "penetrance": "moderate",
            "penetrance_score": 0.4,
            "evidence": "monogenic",
            "note": "脑源性神经营养因子;下丘脑VMH/DMH表达;TrkB信号;调控能量消耗和食欲",
            "phenotype_assoc": "BDNF相关肥胖综合征;常显(外显率可变);认知/行为异常、肥胖、发育迟缓",
            "key_domains": "前体肽(aa1-247);成熟BDNF(aa129-247);半胱氨酸结结构",
            "clingen_validity": "Strong"
        },
        {
            "gene": "NTRK2",
            "tier": "mendelian_mod",
            "contribution_score": 0.7,
            "penetrance": "moderate",
            "penetrance_score": 0.4,
            "evidence": "monogenic",
            "note": "TrkB受体;BDNF受体;下丘脑表达;LOF导致严重肥胖和发育迟缓",
            "phenotype_assoc": "NTRK2相关肥胖综合征;常显/新发;食欲亢进、发育迟缓、肌张力低",
            "key_domains": "富亮氨酸重复域;免疫球蛋白样域;激酶域;BDNF结合域",
            "clingen_validity": "Strong"
        },
        {
            "gene": "SH2B1",
            "tier": "mendelian_mod",
            "contribution_score": 0.7,
            "penetrance": "moderate",
            "penetrance_score": 0.4,
            "evidence": "monogenic",
            "note": "SH2B衔接蛋白1;胰岛素/瘦素受体信号;JAK2调控;缺失导致肥胖和胰岛素抵抗",
            "phenotype_assoc": "SH2B1缺失综合征/单基因肥胖;常显;肥胖、胰岛素抵抗、行为异常",
            "key_domains": "PH域(aa1-100);SH2域(aa350-450);dimerization domain;胰岛素受体相互作用域",
            "clingen_validity": "Moderate"
        },
        {
            "gene": "PCSK1",
            "tier": "mendelian_mod",
            "contribution_score": 0.7,
            "penetrance": "moderate",
            "penetrance_score": 0.4,
            "evidence": "monogenic",
            "note": "前蛋白转化酶1/3;神经内分泌肽加工;POMC/胰岛素原加工;LOF导致肥胖和肠病",
            "phenotype_assoc": "PCSK1缺乏;常隐/复合杂合;早发肥胖、慢性腹泻、内分泌异常",
            "key_domains": "信号肽(aa1-26);前肽(aa27-83);催化域(aa235-596);P域(aa597-753)",
            "clingen_validity": "Strong"
        },
        {
            "gene": "MRAP2",
            "tier": "mendelian_mod",
            "contribution_score": 0.7,
            "penetrance": "moderate",
            "penetrance_score": 0.4,
            "evidence": "monogenic",
            "note": "黑素皮质素受体辅助蛋白2;MC4R辅助蛋白;调控MC4R信号强度",
            "phenotype_assoc": "MRAP2相关肥胖;常显;严重肥胖、肾上腺功能不全",
            "key_domains": "单次跨膜结构域;C端尾巴;MC4R相互作用域",
            "clingen_validity": "Moderate"
        },
        {
            "gene": "CEP19",
            "tier": "mendelian_mod",
            "contribution_score": 0.7,
            "penetrance": "moderate",
            "penetrance_score": 0.4,
            "evidence": "monogenic",
            "note": "中心体蛋白19;调控纤毛/中心体;与mTORC1/胰岛素信号相关",
            "phenotype_assoc": "CEP19相关肥胖;常隐/复合杂合;早发肥胖、轻度智力障碍",
            "key_domains": "CEP19结构域;中心体定位序列",
            "clingen_validity": "Limited"
        },
        {
            "gene": "FTO",
            "tier": "strong_gwas",
            "contribution_score": 0.3,
            "penetrance": "low",
            "penetrance_score": 0.1,
            "evidence": "gwas",
            "note": "BMI/T2D最强GWAS位点;IRX3/IRX5调控;脂肪前体细胞产热",
            "phenotype_assoc": "BMI GWAS最强信号;通过IRX3/IRX5影响脂肪产热;T2D风险部分通过BMI中介",
            "key_domains": "AlkB同源域;2-OG-Fe(II)双加氧酶域;核定位信号",
            "clingen_validity": "Moderate (risk factor)"
        },
        {
            "gene": "TMEM18",
            "tier": "strong_gwas",
            "contribution_score": 0.3,
            "penetrance": "low",
            "penetrance_score": 0.1,
            "evidence": "gwas",
            "note": "跨膜蛋白18;下丘脑高表达;神经元发育/轴突导向;BMI经典位点",
            "phenotype_assoc": "BMI GWAS经典位点;神经元发育相关;儿童期肥胖风险",
            "key_domains": "3次跨膜结构域;N端胞质域;C端胞质域",
            "clingen_validity": ""
        },
        {
            "gene": "GNPDA2",
            "tier": "strong_gwas",
            "contribution_score": 0.3,
            "penetrance": "low",
            "penetrance_score": 0.1,
            "evidence": "gwas",
            "note": "葡糖胺-6-磷酸脱氨酶2;下丘脑表达;能量稳态调控",
            "phenotype_assoc": "BMI GWAS经典位点;hypothalamic nutrient sensing;与胰岛素分泌相关",
            "key_domains": "NADP结合域;催化域",
            "clingen_validity": ""
        },
        {
            "gene": "NEGR1",
            "tier": "strong_gwas",
            "contribution_score": 0.3,
            "penetrance": "low",
            "penetrance_score": 0.1,
            "evidence": "gwas",
            "note": "神经元生长调节因子1;神经细胞黏附;下丘脑表达;BMI经典位点",
            "phenotype_assoc": "BMI GWAS经典位点;神经突生长调控;与食欲/体重调节相关",
            "key_domains": "Ig-like C2型域;Ig-like V型域;FN3域",
            "clingen_validity": ""
        },
        {
            "gene": "MTCH2",
            "tier": "strong_gwas",
            "contribution_score": 0.3,
            "penetrance": "low",
            "penetrance_score": 0.1,
            "evidence": "gwas",
            "note": "线粒体载体同源物2;线粒体功能/凋亡;下丘脑表达;BMI经典位点",
            "phenotype_assoc": "BMI GWAS经典位点;线粒体代谢;脂肪凋亡调控",
            "key_domains": "线粒体载体结构域;6次跨膜",
            "clingen_validity": ""
        },
        {
            "gene": "KLF9",
            "tier": "strong_gwas",
            "contribution_score": 0.3,
            "penetrance": "low",
            "penetrance_score": 0.1,
            "evidence": "gwas",
            "note": "Kruppel样因子9;转录因子;东亚BMI GWAS位点;糖皮质激素响应",
            "phenotype_assoc": "BMI GWAS;东亚人群显著;转录调控;代谢相关",
            "key_domains": "锌指DNA结合域(aa250-330);核定位信号",
            "clingen_validity": ""
        },
        {
            "gene": "LYPLAL1",
            "tier": "gwas",
            "contribution_score": 0.18,
            "penetrance": "low",
            "penetrance_score": 0.1,
            "evidence": "gwas",
            "note": "溶血磷脂酶样1;脂肪代谢;BMI GWAS位点",
            "phenotype_assoc": "BMI GWAS位点;脂肪代谢/溶血磷脂酸调控",
            "key_domains": "patatin-like磷脂酶域;催化三联体",
            "clingen_validity": ""
        },
        {
            "gene": "GPRC5B",
            "tier": "gwas",
            "contribution_score": 0.18,
            "penetrance": "low",
            "penetrance_score": 0.1,
            "evidence": "gwas",
            "note": "G蛋白偶联受体C类5B;脑/脂肪表达;代谢调控",
            "phenotype_assoc": "BMI GWAS位点;GPCR信号;脂肪分化相关",
            "key_domains": "7次跨膜GPCR域;VFT-like胞外域",
            "clingen_validity": ""
        },
        {
            "gene": "FAIM2",
            "tier": "gwas",
            "contribution_score": 0.18,
            "penetrance": "low",
            "penetrance_score": 0.1,
            "evidence": "gwas",
            "note": "Fas凋亡抑制分子2;神经元凋亡调控;BMI经典位点",
            "phenotype_assoc": "BMI GWAS经典位点;神经凋亡调控;下丘脑相关",
            "key_domains": "BIR-like域;Fas结合域",
            "clingen_validity": ""
        },
        {
            "gene": "NRXN3",
            "tier": "gwas",
            "contribution_score": 0.18,
            "penetrance": "low",
            "penetrance_score": 0.1,
            "evidence": "gwas",
            "note": "神经外素3;突触黏附分子;下丘脑表达;BMI GWAS位点",
            "phenotype_assoc": "BMI GWAS位点;突触形成/神经回路;食欲调控",
            "key_domains": "Laminin G域;EGF-like域;PDZ结合域",
            "clingen_validity": ""
        },
        {
            "gene": "SLC39A8",
            "tier": "gwas",
            "contribution_score": 0.18,
            "penetrance": "low",
            "penetrance_score": 0.1,
            "evidence": "gwas",
            "note": "锌转运蛋白ZIP8;全身锌稳态;免疫/代谢多效性",
            "phenotype_assoc": "BMI GWAS位点;锌稳态;免疫与代谢多效性",
            "key_domains": "8次跨膜锌转运域;金属离子结合域",
            "clingen_validity": ""
        },
        {
            "gene": "MAP2K5",
            "tier": "gwas",
            "contribution_score": 0.18,
            "penetrance": "low",
            "penetrance_score": 0.1,
            "evidence": "gwas",
            "note": "MAP激酶激酶5;ERK5通路;脂肪/肌肉代谢;BMI GWAS位点",
            "phenotype_assoc": "BMI GWAS位点;ERK5信号;脂肪代谢",
            "key_domains": "激酶域(aa370-650);PB1域;N端调控域",
            "clingen_validity": ""
        },
        {
            "gene": "QPCTL",
            "tier": "gwas",
            "contribution_score": 0.18,
            "penetrance": "low",
            "penetrance_score": 0.1,
            "evidence": "gwas",
            "note": "谷氨酰肽环转移酶样;谷氨酸环化;神经肽加工;BMI GWAS位点",
            "phenotype_assoc": "BMI GWAS位点;神经肽谷氨酰环化;食欲调控",
            "key_domains": "谷氨酰环转移酶催化域;分泌信号肽",
            "clingen_validity": ""
        },
        {
            "gene": "LRRN6C",
            "tier": "gwas",
            "contribution_score": 0.18,
            "penetrance": "low",
            "penetrance_score": 0.1,
            "evidence": "gwas",
            "note": "富亮氨酸重复神经元6C;神经细胞黏附;下丘脑表达",
            "phenotype_assoc": "BMI GWAS位点;神经元黏附/突触;体重调控",
            "key_domains": "富亮氨酸重复域;Ig-like域;FN3域",
            "clingen_validity": ""
        },
        {
            "gene": "RPL27A",
            "tier": "gwas",
            "contribution_score": 0.15,
            "penetrance": "low",
            "penetrance_score": 0.1,
            "evidence": "gwas",
            "note": "核糖体蛋白L27a;核糖体大亚基;翻译调控",
            "phenotype_assoc": "BMI GWAS位点;核糖体功能;可能通过翻译调控影响体重",
            "key_domains": "核糖体蛋白L27a折叠域",
            "clingen_validity": ""
        },
        {
            "gene": "TNNI3K",
            "tier": "gwas",
            "contribution_score": 0.15,
            "penetrance": "low",
            "penetrance_score": 0.1,
            "evidence": "gwas",
            "note": "心肌肌钙蛋白I相互作用激酶;心脏/代谢;BMI GWAS位点",
            "phenotype_assoc": "BMI GWAS位点;激酶信号;心脏代谢关联",
            "key_domains": "激酶域;N端调节域",
            "clingen_validity": ""
        },
        {
            "gene": "AIF1",
            "tier": "gwas",
            "contribution_score": 0.15,
            "penetrance": "low",
            "penetrance_score": 0.1,
            "evidence": "gwas",
            "note": "同种异体移植炎症因子1;小胶质细胞/巨噬细胞;炎症;代谢",
            "phenotype_assoc": "BMI GWAS位点;炎症调控;免疫-代谢轴",
            "key_domains": "EF-hand Ca2+结合域;免疫球蛋白样域",
            "clingen_validity": ""
        },
        {
            "gene": "FLJ35779",
            "tier": "gwas",
            "contribution_score": 0.15,
            "penetrance": "low",
            "penetrance_score": 0.1,
            "evidence": "gwas",
            "note": "长链非编码RNA/假基因注释;BMI GWAS位点;功能待明确",
            "phenotype_assoc": "BMI GWAS位点;非编码RNA;调控机制待明确",
            "key_domains": "",
            "clingen_validity": ""
        },
        {
            "gene": "CDKAL1",
            "tier": "gwas",
            "contribution_score": 0.18,
            "penetrance": "low",
            "penetrance_score": 0.1,
            "evidence": "gwas",
            "note": "CDK5调节亚基相关蛋白1样1;胰岛β细胞功能/胰岛素分泌;东亚BMI特异性位点",
            "phenotype_assoc": "东亚BMI GWAS位点;T2D与出生体重多效性;胰岛素分泌调控",
            "key_domains": "tRNA甲基转移酶样域;TRAM域",
            "clingen_validity": ""
        },
        {
            "gene": "SEC16B",
            "tier": "gwas",
            "contribution_score": 0.18,
            "penetrance": "low",
            "penetrance_score": 0.1,
            "evidence": "gwas",
            "note": "SEC16同源物B;COPII囊泡形成;ER-高尔基体转运;经典BMI位点",
            "phenotype_assoc": "BMI GWAS经典位点;囊泡转运;内分泌分泌调控",
            "key_domains": "SEC16保守域;COPII结合域",
            "clingen_validity": ""
        },
        {
            "gene": "ADCY3",
            "tier": "gwas",
            "contribution_score": 0.18,
            "penetrance": "low",
            "penetrance_score": 0.1,
            "evidence": "gwas",
            "note": "腺苷酸环化酶3;cAMP信号;下丘脑/脂肪表达;东亚BMI GWAS位点(ADCY3/RBJ区)",
            "phenotype_assoc": "BMI GWAS位点;cAMP信号;能量稳态调控",
            "key_domains": "12次跨膜域;催化域;Gs蛋白结合域",
            "clingen_validity": ""
        },
        {
            "gene": "GIPR",
            "tier": "gwas",
            "contribution_score": 0.18,
            "penetrance": "low",
            "penetrance_score": 0.1,
            "evidence": "gwas",
            "note": "葡萄糖依赖性促胰岛素多肽受体;肠促胰岛素信号;胰腺/脂肪表达",
            "phenotype_assoc": "BMI/肥胖GWAS位点;GIP信号;餐后胰岛素分泌",
            "key_domains": "7次跨膜GPCR域;GIP结合域",
            "clingen_validity": ""
        },
        {
            "gene": "GP2",
            "tier": "gwas",
            "contribution_score": 0.15,
            "penetrance": "low",
            "penetrance_score": 0.1,
            "evidence": "gwas",
            "note": "糖蛋白2;胰腺导管细胞;胰液分泌;东亚BMI GWAS位点",
            "phenotype_assoc": "东亚BMI GWAS位点;胰腺外分泌;代谢调控",
            "key_domains": "C型凝集素样域;跨膜域",
            "clingen_validity": ""
        },
        {
            "gene": "BSN",
            "tier": "gwas",
            "contribution_score": 0.15,
            "penetrance": "low",
            "penetrance_score": 0.1,
            "evidence": "gwas",
            "note": "Bassoon;突触前细胞基质蛋白;神经回路/食欲调控;2025跨族裔肥胖新基因",
            "phenotype_assoc": "跨族裔罕见变异肥胖关联;突触功能;食欲调控",
            "key_domains": "PDZ域;锌指域;coiled-coil域",
            "clingen_validity": ""
        },
        {
            "gene": "YLPM1",
            "tier": "gwas",
            "contribution_score": 0.15,
            "penetrance": "low",
            "penetrance_score": 0.1,
            "evidence": "gwas",
            "note": "YLP motif containing 1;RNA加工/剪接;2025跨族裔肥胖新基因",
            "phenotype_assoc": "跨族裔罕见变异肥胖关联;RNA代谢调控",
            "key_domains": "YLP motif;富含脯氨酸域",
            "clingen_validity": ""
        },
        {
            "gene": "GIGYF1",
            "tier": "gwas",
            "contribution_score": 0.15,
            "penetrance": "low",
            "penetrance_score": 0.1,
            "evidence": "gwas",
            "note": "GRB10 interacting GYF protein 1;胰岛素/IGF信号;2025跨族裔肥胖新基因",
            "phenotype_assoc": "跨族裔罕见变异肥胖关联;胰岛素信号调控",
            "key_domains": "GYF域;PH域;PTB域",
            "clingen_validity": ""
        },
        {
            "gene": "SLC5A3",
            "tier": "gwas",
            "contribution_score": 0.15,
            "penetrance": "low",
            "penetrance_score": 0.1,
            "evidence": "gwas",
            "note": "溶质载体家族5成员3;肌醇转运;脑/外周表达;2025跨族裔肥胖新基因",
            "phenotype_assoc": "跨族裔罕见变异肥胖关联;肌醇稳态;神经代谢",
            "key_domains": "钠/肌醇共转运域",
            "clingen_validity": ""
        },
        {
            "gene": "GRM7",
            "tier": "gwas",
            "contribution_score": 0.15,
            "penetrance": "low",
            "penetrance_score": 0.1,
            "evidence": "gwas",
            "note": "谷氨酸代谢型受体7;GABA/谷氨酸信号;下丘脑/神经回路;2025跨族裔肥胖新基因",
            "phenotype_assoc": "跨族裔罕见变异肥胖关联;神经递质信号;食欲调控",
            "key_domains": "7次跨膜GPCR域;谷氨酸结合域",
            "clingen_validity": ""
        }
    ],
    "known_pathogenic_variants": [
        {
            "rsid": "rs772393451",
            "chrom": "chr18",
            "pos": 60371601,
            "ref": "T",
            "alt": "A",
            "gene": "MC4R",
            "effect_allele": "A",
            "other_allele": "T",
            "beta": None,
            "or_value": None,
            "variant_class": "known_pathogenic",
            "contribution_score": 0.9,
            "confidence": "high",
            "note": "MC4R p.Leu250Gln (c.749T>A); 单基因肥胖综合征; 常显; 食欲亢进; ClinVar Likely Pathogenic"
        },
        {
            "rsid": "rs13447324",
            "chrom": "chr18",
            "pos": 60347089,
            "ref": "C",
            "alt": "A",
            "gene": "MC4R",
            "effect_allele": "A",
            "other_allele": "C",
            "beta": None,
            "or_value": None,
            "variant_class": "known_pathogenic",
            "contribution_score": 0.95,
            "confidence": "high",
            "note": "MC4R p.Tyr35Ter (c.105C>A); 无义变异; 单基因肥胖; 功能缺失; ClinVar Pathogenic"
        },
        {
            "rsid": "rs771139087",
            "chrom": "chr7",
            "pos": 128254657,
            "ref": "G",
            "alt": "-",
            "gene": "LEP",
            "effect_allele": "-",
            "other_allele": "G",
            "beta": None,
            "or_value": None,
            "variant_class": "known_pathogenic",
            "contribution_score": 1.0,
            "confidence": "high",
            "note": "LEP c.398del (p.Gly133fs); 1bp缺失; 先天性瘦素缺乏; 常隐; 严重早发肥胖; ClinVar Pathogenic"
        }
    ],
    "gwas_lead_snps": [
        {
            "rsid": "rs1555543",
            "chrom": "chr1",
            "pos": 96479241,
            "gene": "KLF9",
            "effect_allele": "A",
            "other_allele": "C",
            "beta": -0.0238,
            "or_value": None,
            "eaf_eur": 0.413,
            "eaf_eas": 0.1181,
            "tier": "gwas",
            "variant_class": "gwas_lead",
            "contribution_score": 0.2,
            "confidence": "high",
            "note": "KLF9;东亚BMI显著位点;转录因子;beta=-0.0238 per A allele"
        },
        {
            "rsid": "rs10838738",
            "chrom": "chr11",
            "pos": 47641497,
            "gene": "MTCH2",
            "effect_allele": "A",
            "other_allele": "G",
            "beta": -0.0241,
            "or_value": None,
            "eaf_eur": 0.6452,
            "eaf_eas": 0.6902,
            "tier": "gwas",
            "variant_class": "gwas_lead",
            "contribution_score": 0.2,
            "confidence": "high",
            "note": "MTCH2;线粒体载体;BMI经典位点;beta=-0.0241 per A allele"
        },
        {
            "rsid": "rs2568958",
            "chrom": "chr1",
            "pos": 72299433,
            "gene": "NEGR1",
            "effect_allele": "G",
            "other_allele": "A",
            "beta": -0.0326,
            "or_value": None,
            "eaf_eur": 0.3777,
            "eaf_eas": 0.07698,
            "tier": "gwas",
            "variant_class": "gwas_lead",
            "contribution_score": 0.25,
            "confidence": "high",
            "note": "NEGR1;神经细胞黏附;BMI经典位点;beta=-0.0326 per G allele"
        },
        {
            "rsid": "rs7498665",
            "chrom": "chr16",
            "pos": 28871920,
            "gene": "SH2B1",
            "effect_allele": "G",
            "other_allele": "A",
            "beta": 0.0307,
            "or_value": None,
            "eaf_eur": 0.395,
            "eaf_eas": 0.1243,
            "tier": "gwas",
            "variant_class": "gwas_lead",
            "contribution_score": 0.25,
            "confidence": "high",
            "note": "SH2B1;胰岛素/瘦素信号衔接蛋白;BMI经典位点;beta=0.0307 per G allele"
        },
        {
            "rsid": "rs12444979",
            "chrom": "chr16",
            "pos": 19922278,
            "gene": "GPRC5B",
            "effect_allele": "T",
            "other_allele": "C",
            "beta": -0.0396,
            "or_value": None,
            "eaf_eur": 0.1406,
            "eaf_eas": 0.001354,
            "tier": "gwas",
            "variant_class": "gwas_lead",
            "contribution_score": 0.28,
            "confidence": "high",
            "note": "GPRC5B;G蛋白偶联受体;BMI经典位点;beta=-0.0396 per T allele"
        },
        {
            "rsid": "rs9914578",
            "chrom": "chr17",
            "pos": 2101842,
            "gene": "CEP19",
            "effect_allele": "G",
            "other_allele": "C",
            "beta": 0.0201,
            "or_value": None,
            "eaf_eur": 0.2031,
            "eaf_eas": 0.2117,
            "tier": "gwas",
            "variant_class": "gwas_lead",
            "contribution_score": 0.18,
            "confidence": "high",
            "note": "CEP19;中心体蛋白;单基因肥胖候选;beta=0.0201 per G allele"
        },
        {
            "rsid": "rs7132908",
            "chrom": "chr12",
            "pos": 49869365,
            "gene": "FAIM2",
            "effect_allele": "A",
            "other_allele": "G",
            "beta": 0.0341,
            "or_value": None,
            "eaf_eur": 0.3935,
            "eaf_eas": 0.2294,
            "tier": "gwas",
            "variant_class": "gwas_lead",
            "contribution_score": 0.25,
            "confidence": "high",
            "note": "FAIM2;Fas凋亡抑制分子;BMI经典位点;beta=0.0341 per A allele"
        },
        {
            "rsid": "rs10146997",
            "chrom": "chr14",
            "pos": 79478819,
            "gene": "NRXN3",
            "effect_allele": "G",
            "other_allele": "A",
            "beta": 0.025,
            "or_value": None,
            "eaf_eur": 0.2198,
            "eaf_eas": 0.000774,
            "tier": "gwas",
            "variant_class": "gwas_lead",
            "contribution_score": 0.2,
            "confidence": "high",
            "note": "NRXN3;神经外素3;突触黏附;BMI位点;beta=0.0250 per G allele"
        },
        {
            "rsid": "rs2112347",
            "chrom": "chr5",
            "pos": 75719417,
            "gene": "FLJ35779",
            "effect_allele": "G",
            "other_allele": "T",
            "beta": -0.0261,
            "or_value": None,
            "eaf_eur": 0.3629,
            "eaf_eas": 0.5537,
            "tier": "gwas",
            "variant_class": "gwas_lead",
            "contribution_score": 0.2,
            "confidence": "high",
            "note": "FLJ35779;BMI GWAS位点;非编码;beta=-0.0261 per G allele"
        },
        {
            "rsid": "rs2241423",
            "chrom": "chr15",
            "pos": 67794500,
            "gene": "MAP2K5",
            "effect_allele": "A",
            "other_allele": "G",
            "beta": -0.031,
            "or_value": None,
            "eaf_eur": 0.2194,
            "eaf_eas": 0.5941,
            "tier": "gwas",
            "variant_class": "gwas_lead",
            "contribution_score": 0.25,
            "confidence": "high",
            "note": "MAP2K5;ERK5激酶;BMI位点;beta=-0.0310 per A allele"
        },
        {
            "rsid": "rs2287019",
            "chrom": "chr19",
            "pos": 45698914,
            "gene": "QPCTL",
            "effect_allele": "C",
            "other_allele": "T",
            "beta": 0.036,
            "or_value": None,
            "eaf_eur": 0.8074,
            "eaf_eas": 0.8156,
            "tier": "gwas",
            "variant_class": "gwas_lead",
            "contribution_score": 0.28,
            "confidence": "high",
            "note": "QPCTL;谷氨酰肽环转移酶;BMI位点;beta=0.0360 per C allele"
        },
        {
            "rsid": "rs10968576",
            "chrom": "chr9",
            "pos": 28414341,
            "gene": "LRRN6C",
            "effect_allele": "G",
            "other_allele": "A",
            "beta": 0.0249,
            "or_value": None,
            "eaf_eur": 0.3158,
            "eaf_eas": 0.1772,
            "tier": "gwas",
            "variant_class": "gwas_lead",
            "contribution_score": 0.2,
            "confidence": "high",
            "note": "LRRN6C;富亮氨酸重复神经元蛋白;BMI位点;beta=0.0249 per G allele"
        },
        {
            "rsid": "rs1514175",
            "chrom": "chr1",
            "pos": 74525960,
            "gene": "TNNI3K",
            "effect_allele": "A",
            "other_allele": "G",
            "beta": 0.023,
            "or_value": None,
            "eaf_eur": 0.4279,
            "eaf_eas": 0.7564,
            "tier": "gwas",
            "variant_class": "gwas_lead",
            "contribution_score": 0.2,
            "confidence": "high",
            "note": "TNNI3K;心肌肌钙蛋白I相互作用激酶;BMI位点;beta=0.0230 per A allele"
        },
        {
            "rsid": "rs4929949",
            "chrom": "chr11",
            "pos": 8583046,
            "gene": "RPL27A",
            "effect_allele": "C",
            "other_allele": "T",
            "beta": 0.0173,
            "or_value": None,
            "eaf_eur": 0.5026,
            "eaf_eas": 0.4313,
            "tier": "gwas",
            "variant_class": "gwas_lead",
            "contribution_score": 0.15,
            "confidence": "high",
            "note": "RPL27A;核糖体蛋白;BMI位点;beta=0.0173 per C allele"
        },
        {
            "rsid": "rs13107325",
            "chrom": "chr4",
            "pos": 102267552,
            "gene": "SLC39A8",
            "effect_allele": "C",
            "other_allele": "T",
            "beta": -0.0477,
            "or_value": None,
            "eaf_eur": 0.9296,
            "eaf_eas": 0.9996,
            "tier": "gwas",
            "variant_class": "gwas_lead",
            "contribution_score": 0.3,
            "confidence": "high",
            "note": "SLC39A8;锌转运蛋白ZIP8;BMI位点;免疫-代谢多效性;beta=-0.0477 per C allele"
        },
        {
            "rsid": "rs2206734",
            "chrom": "chr6",
            "pos": 20694653,
            "gene": "CDKAL1",
            "effect_allele": "T",
            "other_allele": "C",
            "beta": 0.015,
            "or_value": None,
            "eaf_eur": 0.15,
            "eaf_eas": 0.12,
            "tier": "gwas",
            "variant_class": "gwas_lead",
            "contribution_score": 0.18,
            "confidence": "high",
            "note": "CDKAL1;Wen 2012东亚BMI新位点;胰岛β细胞/胰岛素分泌;EAS中T allele升高BMI;beta≈+0.015 per T allele"
        },
        {
            "rsid": "rs261967",
            "chrom": "chr5",
            "pos": 96514546,
            "gene": "PCSK1",
            "effect_allele": "A",
            "other_allele": "C",
            "beta": 0.015,
            "or_value": None,
            "eaf_eur": 0.5417,
            "eaf_eas": 0.35,
            "tier": "gwas",
            "variant_class": "gwas_lead",
            "contribution_score": 0.18,
            "confidence": "high",
            "note": "PCSK1;Wen 2012东亚BMI新位点;前蛋白转化酶;EAS中A allele升高BMI;beta≈+0.015 per A allele"
        },
        {
            "rsid": "rs2605100",
            "chrom": "chr1",
            "pos": 219470882,
            "gene": "SEC16B",
            "effect_allele": "A",
            "other_allele": "G",
            "beta": 0.009,
            "or_value": None,
            "eaf_eur": 0.3083,
            "eaf_eas": 0.25,
            "tier": "gwas",
            "variant_class": "gwas_lead",
            "contribution_score": 0.18,
            "confidence": "high",
            "note": "SEC16B;COPII囊泡转运;经典BMI位点;GIANT EUR beta=+0.009 per A allele"
        }
    ],
    "prs_variants": [],
    "prs_variants_high": [
        {
            "rsid": "rs9939609",
            "chrom": "chr16",
            "pos": 53786615,
            "gene": "FTO",
            "effect_allele": "A",
            "other_allele": "T",
            "beta": 0.0776,
            "or_value": None,
            "eaf_eur": 0.4079,
            "eaf_eas": 0.1548,
            "tier": "prs_high",
            "variant_class": "prs_high",
            "contribution_score": 0.45,
            "confidence": "high",
            "note": "FTO;BMI最强GWAS位点;IRX3/IRX5调控;beta=0.0776 per A allele"
        },
        {
            "rsid": "rs17782313",
            "chrom": "chr18",
            "pos": 60183864,
            "gene": "MC4R",
            "effect_allele": "C",
            "other_allele": "T",
            "beta": 0.0566,
            "or_value": None,
            "eaf_eur": 0.2348,
            "eaf_eas": 0.1833,
            "tier": "prs_high",
            "variant_class": "prs_high",
            "contribution_score": 0.4,
            "confidence": "high",
            "note": "MC4R;BMI经典位点;MC4R下游调控;beta=0.0566 per C allele"
        },
        {
            "rsid": "rs10767664",
            "chrom": "chr11",
            "pos": 27704439,
            "gene": "BDNF",
            "effect_allele": "A",
            "other_allele": "T",
            "beta": 0.043,
            "or_value": None,
            "eaf_eur": 0.7892,
            "eaf_eas": 0.5301,
            "tier": "prs_high",
            "variant_class": "prs_high",
            "contribution_score": 0.35,
            "confidence": "high",
            "note": "BDNF;脑源性神经营养因子;下丘脑;beta=0.0430 per A allele"
        },
        {
            "rsid": "rs2867125",
            "chrom": "chr2",
            "pos": 622827,
            "gene": "TMEM18",
            "effect_allele": "T",
            "other_allele": "C",
            "beta": -0.0592,
            "or_value": None,
            "eaf_eur": 0.1721,
            "eaf_eas": 0.08073,
            "tier": "prs_high",
            "variant_class": "prs_high",
            "contribution_score": 0.35,
            "confidence": "high",
            "note": "TMEM18;BMI经典位点;TMEM18 T allele降低BMI(beta=-0.0592);C allele增加"
        },
        {
            "rsid": "rs10938397",
            "chrom": "chr4",
            "pos": 45180510,
            "gene": "GNPDA2",
            "effect_allele": "A",
            "other_allele": "G",
            "beta": -0.0402,
            "or_value": None,
            "eaf_eur": 0.5691,
            "eaf_eas": 0.7175,
            "tier": "prs_high",
            "variant_class": "prs_high",
            "contribution_score": 0.3,
            "confidence": "high",
            "note": "GNPDA2;BMI经典位点;beta=-0.0402 per A allele"
        }
    ],
    "regulatory_regions": [
        {
            "chrom": "chr16",
            "start": 53750000,
            "end": 53810000,
            "gene": "FTO",
            "type": "gwas_locus",
            "contribution_score": 0.45,
            "note": "FTO/IRX3/IRX5 BMI最强GWAS区;脂肪产热调控"
        },
        {
            "chrom": "chr2",
            "start": 600000,
            "end": 650000,
            "gene": "TMEM18",
            "type": "gwas_locus",
            "contribution_score": 0.35,
            "note": "TMEM18 BMI经典位点;下丘脑神经元发育"
        },
        {
            "chrom": "chr18",
            "start": 60150000,
            "end": 60200000,
            "gene": "MC4R",
            "type": "gwas_locus",
            "contribution_score": 0.4,
            "note": "MC4R BMI经典位点/单基因肥胖区"
        },
        {
            "chrom": "chr11",
            "start": 27650000,
            "end": 27750000,
            "gene": "BDNF",
            "type": "gwas_locus",
            "contribution_score": 0.3,
            "note": "BDNF BMI GWAS区;下丘脑能量稳态"
        }
    ],
    "key_regions": {
        "MC4R": {
            "note": "黑素皮质素4受体;7次跨膜GPCR;下丘脑食欲/能量稳态核心",
            "regions": [
                {
                    "name": "Domain",
                    "residues": "1-40",
                    "note": "N端胞外域"
                },
                {
                    "name": "Domain",
                    "residues": "42-130",
                    "note": "配体结合域"
                },
                {
                    "name": "Domain",
                    "residues": "220-333",
                    "note": "G蛋白偶联域"
                }
            ],
            "critical_residues": []
        },
        "LEP": {
            "note": "瘦素;脂肪信号;JAK-STAT通路",
            "regions": [
                {
                    "name": "Domain",
                    "residues": "1-21",
                    "note": "信号肽"
                },
                {
                    "name": "Domain",
                    "residues": "22-167",
                    "note": "核心瘦素"
                }
            ],
            "critical_residues": []
        },
        "POMC": {
            "note": "阿黑皮素原;下丘脑/垂体;裂解为α-MSH/ACTH/β-内啡肽",
            "regions": [
                {
                    "name": "Domain",
                    "residues": "1-26",
                    "note": "信号肽"
                },
                {
                    "name": "Domain",
                    "residues": "138-176",
                    "note": "ACTH"
                },
                {
                    "name": "Domain",
                    "residues": "138-150",
                    "note": "α-MSH"
                }
            ],
            "critical_residues": []
        },
        "BDNF": {
            "note": "脑源性神经营养因子;TrkB受体;下丘脑VMH/DMH;能量消耗调控",
            "regions": [
                {
                    "name": "Domain",
                    "residues": "1-247",
                    "note": "前体肽"
                },
                {
                    "name": "Domain",
                    "residues": "129-247",
                    "note": "成熟BDNF"
                }
            ],
            "critical_residues": []
        },
        "FTO": {
            "note": "FTO;IRX3/IRX5调控;脂肪产热;BMI最强GWAS位点",
            "regions": [
                {
                    "name": "Domain",
                    "residues": "1-500",
                    "note": "AlkB双加氧酶域"
                }
            ],
            "critical_residues": []
        },
        "SH2B1": {
            "note": "SH2B1;胰岛素/瘦素受体信号衔接蛋白",
            "regions": [
                {
                    "name": "Domain",
                    "residues": "1-100",
                    "note": "PH域"
                },
                {
                    "name": "Domain",
                    "residues": "350-450",
                    "note": "SH2域"
                }
            ],
            "critical_residues": []
        },
        "PCSK1": {
            "note": "PCSK1;神经内分泌肽前体加工酶;POMC/胰岛素原加工",
            "regions": [
                {
                    "name": "Domain",
                    "residues": "1-26",
                    "note": "信号肽"
                },
                {
                    "name": "Domain",
                    "residues": "235-596",
                    "note": "催化域"
                },
                {
                    "name": "Domain",
                    "residues": "597-753",
                    "note": "P域"
                }
            ],
            "critical_residues": []
        }
    },
    "key_literature": [
        {
            "pmid": "25673413",
            "title": "Genetic studies of body mass index yield new insights for obesity biology",
            "genes": [
                "FTO",
                "MC4R",
                "BDNF",
                "TMEM18",
                "GNPDA2",
                "NEGR1",
                "MTCH2",
                "SH2B1"
            ],
            "note": "Locke et al., Nature 2015; 339,224人EUR BMI GWAS; 97个BMI位点; CNS在肥胖易感性中的核心作用",
            "evidence_type": "gwas"
        },
        {
            "pmid": "22344221",
            "title": "Genome-wide association study identifies six new loci influencing body mass index in East Asians",
            "genes": [
                "FTO",
                "MC4R",
                "BDNF",
                "CDKAL1",
                "PCSK1",
                "SEC16B",
                "ADCY3",
                "GIPR",
                "GP2",
                "KLF9",
                "GPRC5B"
            ],
            "note": "Wen et al., Nature Genetics 2012; 东亚BMI GWAS; 发现CDKAL1/PCSK1/SEC16B等东亚显著位点",
            "evidence_type": "gwas"
        },
        {
            "pmid": "39506477",
            "title": "Discovery of obesity genes through cross-ancestry analysis",
            "genes": [
                "YLPM1",
                "GIGYF1",
                "SLC5A3",
                "GRM7",
                "BSN",
                "MC4R"
            ],
            "note": "Banerjee & Girirajan, Nature Communications 2025; 83.9万跨族裔罕见变异分析; 13个BMI相关基因",
            "evidence_type": "gwas"
        },
        {
            "pmid": "21921914",
            "title": "Melanocortin-4 receptor mutations in obesity",
            "genes": [
                "MC4R",
                "LEP",
                "LEPR",
                "POMC",
                "PCSK1"
            ],
            "note": "Huszar et al./review; 单基因肥胖综合征; MC4R通路核心",
            "evidence_type": "review"
        }
    ]
}
