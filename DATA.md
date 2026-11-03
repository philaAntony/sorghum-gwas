# Data description

## Source
Dryad dataset doi:10.5061/dryad.3j9kd51w9 (public domain): "Integrating SNP, RNASeq, and Phenotypic data from sorghum association and diversity panels for comprehensive genomic analysis" (Subash, Mangal, Shrestha, Xu, Schnable, Yang; University of Nebraska-Lincoln). Raw files are not stored in this repository (`data/` is git-ignored).

## Genotypes
File: `B.SDP_SAP_815_maf001_het10_imputed.recode.vcf.gz`
- 815 lines (Sorghum Association Panel plus Sorghum Diversity Panel), SNPs called from RNA-seq of mature leaf tissue against the BTx623 v5 reference (as stated in the dataset README; not independently verified).
- 346,222 variants: 277,724 SNPs (all biallelic) and 68,498 indels; chromosomes Chr01 to Chr10 only.
- Variants with MAF < 0.01 or > 10% heterozygous calls were removed by the data authors. The VCF header shows Beagle processing; no marker carries an imputation flag or INFO values, so imputed genotypes cannot be identified or filtered by imputation quality.
- No missing genotypes (0 of 282,170,930 calls); 1.07% of calls are heterozygous.
- Markers come from expressed regions, so genome coverage is uneven.

## Phenotypes
File: 2021 Lincoln, Nebraska field experiment (`2.SDP_SAP_f2021_field_phenotype`), high-nitrogen (HN) treatment only.
- Traits: plant height and days to flower, averaged per line over replicates (`metadata/pheno_2021_HN.tsv`, produced by `scripts/01_make_phenotype_table.py`).
- Plant height units: the dataset README says centimetres, but the values (median 0.94) indicate metres. Units are recorded here as metres (inferred).
- Two implausible plant height records (6.0 for PI656121 and 9.67 for PI534063) are treated as missing (rule: height > 5 m). Heights up to 4.5 m were kept; the distribution is revisited in Week 11.
- Other years were not used: files differ in design and scale (for example the 2023 median height is 251 against 95.5 in 2020, and 2023 contains implausible values).

## Matching genotypes to phenotypes
Phenotype IDs contain spaces and varying case (for example "PI 152651"); the VCF uses names without spaces ("PI152651"). IDs are matched after removing spaces, hyphens and underscores and ignoring case, using the PI number, then the conversion ID, then the line name (scripts 00b and 01). In the 2021 file all 815 VCF lines have a matching entry.

## Known issues and open questions
- Replicate structure: high nitrogen has two replicates (4 and 5, 966 plots each) and low nitrogen two (3 and 6, 390 plots each). High nitrogen covers the SAP plus the diversity panel; low nitrogen covers the SAP only. 37 lines have a single usable plant height record, so their means are noisier.
- PI564163 is the check line BTx623, replicated in 190 plots (both treatments); it is kept as a single line. Its mean height of 1.10 supports the metre interpretation.
- Eight VCF sample names are not accession IDs (Border, Day, ElMota, RTAM2566, FeteritaGishesh, SanChiSan, Comb7078, Segaolane); Border is probably a field label and will be removed.
- The reference build (v5) and a matching gene annotation must be confirmed before candidate gene lookup (Week 15).
- The panel mixes SAP and SDP lines, so population structure is expected to be strong.

## Phenotype table summary
- Lines with plant height: 813; with days to flower: 813; with both: 813 (of 815 genotyped lines).

## Quality control (Week 11)
Scripts: `02_genotype_qc.sh`, `03_phenotype_qc.R`, `04_make_gwas_phenotypes.R`. SNP counts: `results/qc_snp_counts.tsv`; summary: `results/qc_summary.txt`.

| Stage | Samples | Variants |
|---|---|---|
| Raw VCF | 815 | 346,222 |
| Biallelic SNPs, kept samples | 812 | 277,724 |
| After MAF >= 0.05 | 812 | 164,212 |
| After LD thinning (r2 < 0.2, 50-SNP window) | 812 | 23,900 |

- **Samples:** 812 lines kept; `Border` (a field label) was removed. `Day` was kept because it has genotypes and plausible phenotypes, but its identity is unclear.
- **MAF 0.05** removes 113,512 SNPs, mostly in the 0.01 to 0.05 range. With 812 lines, rare alleles carry little power; a lower cutoff can be tried as a sensitivity check (`MAF=0.01 bash scripts/02_genotype_qc.sh`).
- **The LD-thinned set (23,900 SNPs) is for PCA and kinship only;** the GWAS uses the full MAF-filtered set (164,212 SNPs).
- **Quality:** Ts/Tv ratio 1.586; no missing genotypes. Per-line heterozygosity has a median of 0.9% and a maximum of 5.6% (PI533938); PI533994, PI576385, SC1442 and PI655978 are also above 3%. These lines were kept and will be checked in the PCA.
- **Plant height** is right-skewed (skewness 2.01; 0.65 after log transform), so the GWAS uses log plant height. 16 lines have |z| > 3, all tall (2.4 to 3.3 m); they were kept and flagged, with a sensitivity analysis planned. Days to flower (skewness 0.38) is used untransformed.
- Height and days to flower are weakly correlated (r = 0.138).
- 37 lines have a single plant height record, so their means are noisier.
- The QQ plot of log plant height follows the line in the centre but has a heavier right tail (the tall lines). A rank-based inverse-normal transform and a rerun without the 16 flagged lines are planned as sensitivity checks.

## Population structure (Week 12)
Scripts: `05_line_groups.py`, `06_population_structure.R`, `07_pc_drivers.R`. Outputs in `results/` and `figures/`. Based on the 23,900 LD-thinned SNPs in 812 lines.
- **Kinship (VanRaden):** diagonal mean 1.94 (range 0.96 to 2.44), as expected for inbred lines; median off-diagonal -0.047, 99th percentile 0.58. The heatmap shows blocks of related lines and near-zero relatedness between blocks.
- **Near-duplicates:** one pair above 0.9 standardised relatedness, PI656084 and PI564163 (the BTx623 check), 0.906. Both kept.
- **PCA:** PC1 5.4%, PC2 3.9%, PC3 3.5%, PC4 1.9%, PC5 1.4% (cumulative 16.1%): weak structure with an elbow after PC3 and a long tail. The PC1-PC2 plot is triangular, consistent with admixture among a few groups (not tested). SAP and non-SAP lines are intermixed; SC-named lines concentrate on the PC1-negative arm. The k-means clusters are exploratory only.
- **PC drivers:** for PC1 to PC6 the top 1% of correlated SNPs (239 per PC) are spread over all ten chromosomes (at most 27% on one chromosome), so these PCs reflect genome-wide structure, not single regions. Not normalised for SNP count per chromosome; chromosome 7 is under-represented for all PCs.
- **Structure versus traits:** PC1 to 5 explain 18.2% of the variance in days to flower (correlations -0.23, -0.30 and 0.18 with PC1 to 3; cluster means 62.9 to 70.3 days) but only 3.3% for log plant height (PC5 r = -0.17).
- **Implications:** flowering time is the trait where uncorrected structure should inflate the GWAS; plant height is the positive-control trait. Planned models: no correction, PC1-3, PC1-5, kinship, kinship + PC1-3, compared by QQ plots and lambda.
- 336 of the 812 lines are in the SAP (the 2021 phenotype file lists 376 SAP lines, so about 40 have no genotypes in this VCF).
