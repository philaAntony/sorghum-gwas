## Week 10: Choosing data and planning the GWAS
- **Did:** switched from the rice panel to sorghum, picked the SAP/SDP data on Dryad, matched genotype and phenotype IDs, built a per-line phenotype table, wrote DATA.md and PLAN.md.
- **Broke:** Dryad refused command-line downloads (403) so I downloaded in the browser; phenotype IDs had spaces and different case from the VCF; the plant height values were in metres although the README said centimetres, and two records were impossible.
- **Learned:** check a dataset's documentation against the data itself; harmonise sample IDs in a script; a repeated check line (BTx623) is useful for sanity checks.

## Week 11: Genotype and phenotype QC
- **Did:** filtered to biallelic SNPs in 812 lines, applied MAF >= 0.05 (277,724 to 164,212 SNPs), thinned for LD (23,900 SNPs) for structure analysis, checked per-line heterozygosity and Ts/Tv, and chose a log transform for plant height.
- **Broke:** nothing major; had to rename chromosomes to numbers for PLINK, and a stray .Rhistory file nearly got committed.
- **Learned:** an outlier is not an error (the tall lines look biological); use the LD-pruned set only for structure, not for the GWAS; keep the thresholds and sample lists in scripts so the QC is reproducible.

## Week 12: Population structure
- **Did:** built a VanRaden kinship matrix and PCA from 23,900 LD-thinned SNPs, checked near-duplicate lines, looked at how the PCs relate to SAP membership, line name types and the two traits.
- **Broke:** I assumed the conversion-ID column marked conversion lines, but 683 of 812 lines had a value, so I dropped it and used the SC name type instead.
- **Learned:** structure here is weak and continuous (PC1 5.4%); flowering time is far more structured than plant height; a low-variance PC can track a trait because of a local region, so check what drives a PC before using it as a covariate.
