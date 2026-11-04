suppressMessages(library(rMVP)); options(rMVP.OutputLog2File = FALSE)
genotype <- attach.big.matrix("data/work/mvp_maf.geno.desc")
K <- attach.big.matrix("data/work/mvp_maf.kin.desc")
map <- read.table("data/work/mvp_maf.geno.map", header = TRUE, stringsAsFactors = FALSE)
ind <- read.table("data/work/mvp_maf.geno.ind", header = FALSE, colClasses = "character")[, 1]
ph <- read.delim("metadata/pheno_gwas.tsv", colClasses = c("character", "numeric", "numeric"))
ph <- ph[match(ind, ph$taxa), ]; stopifnot(identical(ph$taxa, ind))
sc <- read.delim("results/pca_scores.tsv", colClasses = c(taxa = "character"))
pcs <- as.matrix(sc[match(ind, sc$taxa), paste0("PC", 1:3)])
flag <- read.delim("results/phenotype_flagged.tsv", colClasses = c(taxa = "character"))
m <- nrow(map); bonf <- 0.05 / m
rint <- function(x) { r <- rank(x, na.last = "keep"); qnorm((r - 0.5) / sum(!is.na(x))) }
lambda <- function(p) { p <- pmax(p[!is.na(p)], 1e-300); median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1) }

rows <- list()
specs <- list(
  list(tr = "log_plant_height", flagtrait = "plant_height"),
  list(tr = "days_to_flower", flagtrait = "days_to_flower"))
for (s in specs) {
  main <- read.delim(gzfile(sprintf("results/gwas/%s_MLM_K_PC3.tsv.gz", s$tr)))
  lead <- main$snp[which.min(main$p)]
  drop_ids <- flag$taxa[flag$trait == s$flagtrait]
  variants <- list("drop flagged lines" = function(y) { y[ind %in% drop_ids] <- NA; y },
                   "rank-based inverse normal" = rint)
  for (v in names(variants)) {
    y <- variants[[v]](ph[[s$tr]])
    mvp <- MVP(phe = data.frame(taxa = ind, y = y), geno = genotype, map = map, K = K,
               CV.MLM = pcs, method = "MLM", ncpus = 4, file.output = FALSE, verbose = FALSE)
    res <- as.matrix(mvp$mlm.results); p <- res[, ncol(res)]
    rows[[length(rows) + 1]] <- data.frame(trait = s$tr, variant = v, n_lines = sum(!is.na(y)),
      lambda = round(lambda(p), 3), n_bonferroni = sum(p < bonf, na.rm = TRUE),
      top_snp = map[which.min(p), 1], top_p = signif(min(p, na.rm = TRUE), 3),
      p_at_main_lead = signif(p[match(lead, map[, 1])], 3), main_lead = lead,
      spearman_vs_main = round(cor(-log10(pmax(p, 1e-300)), -log10(pmax(main$p, 1e-300)), method = "spearman", use = "complete.obs"), 3))
  }
}
out <- do.call(rbind, rows)
write.table(out, "results/sensitivity_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(out, row.names = FALSE)
