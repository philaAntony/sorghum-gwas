## Week 10: Choosing data and planning the GWAS
- **Did:** switched from the rice panel to sorghum, picked the SAP/SDP data on Dryad, matched genotype and phenotype IDs, built a per-line phenotype table, wrote DATA.md and PLAN.md.
- **Broke:** Dryad refused command-line downloads (403) so I downloaded in the browser; phenotype IDs had spaces and different case from the VCF; the plant height values were in metres although the README said centimetres, and two records were impossible.
- **Learned:** check a dataset's documentation against the data itself; harmonise sample IDs in a script; a repeated check line (BTx623) is useful for sanity checks.

## Week 11: Genotype and phenotype QC
- **Did:** filtered to biallelic SNPs in 812 lines, applied MAF >= 0.05 (277,724 to 164,212 SNPs), thinned for LD (23,900 SNPs) for structure analysis, checked per-line heterozygosity and Ts/Tv, and chose a log transform for plant height.
- **Broke:** nothing major; had to rename chromosomes to numbers for PLINK, and a stray .Rhistory file nearly got committed.
- **Learned:** an outlier is not an error (the tall lines look biological); use the LD-pruned set only for structure, not for the GWAS; keep the thresholds and sample lists in scripts so the QC is reproducible.

## Week 12: Population structure
- **Did:** built a VanRaden kinship matrix and PCA from 23,900 LD-thinned SNPs, checked near-duplicate lines, and related the PCs to SAP membership, name types and the two traits.
- **Broke:** I assumed the conversion-ID column marked conversion lines, but 683 of 812 lines had a value, so I used the SC name type instead.
- **Learned:** structure here is weak and fine-scale (PC1 only 5.4%); flowering time is far more structured than plant height; I guessed a height-linked PC reflected one genomic region, but the check showed it was genome-wide, so test such guesses before acting on them.

## Week 13: GWAS with and without structure correction
- **Did:** converted the MAF-filtered SNPs to rMVP format, ran GLM and mixed-model (MLM) GWAS for log plant height and days to flower, with and without PC1-3, and compared inflation (lambda) across models.
- **Broke:** I predicted plant height would barely be affected by structure, but the naive GLM was inflated (lambda 2.2); I committed before looking at the results.
- **Learned:** kinship, not a few PCs, removes the inflation in this panel; a naive GWAS produces thousands of false hits (19,341 SNPs for flowering time); the structure-driven chr5 flowering hit disappeared once structure was modelled.

## Week 14: Visualising and diagnosing the GWAS
- **Did:** made Manhattan and QQ plots for four models per trait, set significant and suggestive thresholds, listed loci, compared how many hits each model shares, ran sensitivity analyses (dropping outlying lines, inverse-normal transform), and checked LD between the two chromosome 9 lead SNPs.
- **Broke:** nothing major; I had guessed that structure would hardly matter for plant height and was wrong (lambda 2.2 in the naive model).
- **Learned:** kinship, not PCs, fixes the inflation here; a naive GWAS gave 19,341 "hits" for flowering time against 15 in the mixed model; hits just past a threshold are fragile and their effect sizes likely overestimated.
