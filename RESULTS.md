# Results (work in progress)

## GWAS models compared (Week 13)
rMVP; 164,212 SNPs, 812 lines. Bonferroni p = 3.04e-07 (-log10 6.52); 1/m p = 6.09e-06 (-log10 5.22). Full tables are in `results/gwas/` (not tracked); summaries in `results/gwas_lambda.tsv` and `results/gwas_top_hits.tsv`.

| Trait | Model | lambda | SNPs past Bonferroni | SNPs past 1/m | Top SNP |
|---|---|---|---|---|---|
| log plant height | GLM, no PCs | 2.155 | 1230 | 1888 | 9:60843527 |
| log plant height | GLM + PC1-3 | 2.38 | 1447 | 2257 | 9:60843527 |
| log plant height | MLM, kinship | 0.921 | 77 | 128 | 9:60868876 |
| log plant height | MLM, kinship + PC1-3 | 0.919 | 77 | 128 | 9:60868876 |
| days to flower | GLM, no PCs | 8.388 | 19341 | 27498 | 5:55801375 |
| days to flower | GLM + PC1-3 | 2.064 | 153 | 485 | 9:62620720 |
| days to flower | MLM, kinship | 0.966 | 15 | 38 | 9:62620720 |
| days to flower | MLM, kinship + PC1-3 | 0.962 | 15 | 39 | 9:62620720 |

- Mixed models bring lambda close to 1; the naive GLM is strongly inflated, most for days to flower (8.4).
- Adding PC1-3 alone does not fix the inflation; kinship does. With kinship in the model, PCs change almost nothing.
- Plant height: one strong locus on chromosome 9 (about 60.84 to 60.87 Mb) in every model, effect about 0.17 on the log scale.
- Days to flower: the naive GLM's top hit (chr5) is not supported once structure is modelled; the mixed models find chr9:62.62 Mb, about 1.75 Mb from the height peak.
- Positions are on the BTx623 v5 build; no gene annotation yet.
- Caveats: inflation in the GLM is not only structure (a strong locus in LD also raises lambda); lambda slightly below 1 under the mixed models may reflect the kinship matrix including the tested SNP's region.
