# Project 2 plan: GWAS of plant height in sorghum

**Question:** which genomic regions are associated with plant height (and, secondarily, days to flower) in sorghum under high nitrogen, and how much does correcting for population structure change the result?

**Data:** see DATA.md (815 lines, 277,724 SNPs, 2021 Lincoln HN phenotypes).

| Week | Task | Output |
|---|---|---|
| 11 | Genotype QC: keep biallelic SNPs, check MAF and heterozygosity, LD thinning for structure analysis; phenotype QC: distribution, outliers, transformation, remove non-accession lines | filtering script, SNP counts before and after, phenotype summary |
| 12 | Population structure: PCA and kinship matrix; colour the PCA by line type (SC conversion, PI, BTx/RTx) and by PCA clusters | PCA figures, notes on implications |
| 13 | GWAS in GAPIT or rMVP: first without structure correction, then with PCs and kinship (mixed model) | scripts, raw result tables |
| 14 | Manhattan and QQ plots, significance threshold, compare models | figures with captions |
| 15 | Candidate genes near significant SNPs using a v5-matched annotation; check known height genes (Dw1, Dw2, Dw3) as positive controls | candidate gene table with references |
| 16 | README, reproducible workflow, portfolio | repo ready to share |

**Risks to manage:** unit and outlier problems in the phenotype file; SAP plus SDP structure; uneven RNA-seq marker coverage; unknown imputation quality; annotation build mismatch; replicate structure.
