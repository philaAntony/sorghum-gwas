suppressMessages(library(rMVP))
options(rMVP.OutputLog2File = FALSE)
dir.create("results/gwas", recursive = TRUE, showWarnings = FALSE)

# ---- inputs ----
genotype <- attach.big.matrix("data/work/mvp_maf.geno.desc")
K <- attach.big.matrix("data/work/mvp_maf.kin.desc")
map <- read.table("data/work/mvp_maf.geno.map", header = TRUE, stringsAsFactors = FALSE)
ind <- read.table("data/work/mvp_maf.geno.ind", header = FALSE, colClasses = "character")[, 1]
n_snps <- length(readLines("data/work/snps_maf.bim"))
n_ind <- length(readLines("data/work/snps_maf.fam"))
cat("genotype dim:", dim(genotype), "| kinship dim:", dim(K), "| map rows:", nrow(map), "| individuals:", length(ind), "\n")
stopifnot(nrow(map) == n_snps, length(ind) == n_ind)

ph <- read.delim("metadata/pheno_gwas.tsv", colClasses = c("character", "numeric", "numeric"))
ph <- ph[match(ind, ph$taxa), ]
stopifnot(!anyNA(ph$taxa), identical(ph$taxa, ind))      # phenotype order must equal genotype order
sc <- read.delim("results/pca_scores.tsv", colClasses = c(taxa = "character"))
pcs <- as.matrix(sc[match(ind, sc$taxa), paste0("PC", 1:3)])
stopifnot(!anyNA(pcs))

m <- nrow(map); bonf <- 0.05 / m; sugg <- 1 / m
cat(sprintf("SNPs: %d | Bonferroni p = %.2e (-log10 %.2f) | 1/m p = %.2e (-log10 %.2f)\n",
            m, bonf, -log10(bonf), sugg, -log10(sugg)))

lambda <- function(p) {
  p <- pmax(p[!is.na(p)], 1e-300)
  median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
}

summ <- list(); tops <- list(); first <- TRUE
for (tr in c("log_plant_height", "days_to_flower")) {
  y <- data.frame(taxa = ind, y = ph[[tr]])
  for (set in c("none", "PC3")) {
    cv <- if (set == "none") NULL else pcs
    t0 <- proc.time()[3]
    mvp <- MVP(phe = y, geno = genotype, map = map, K = K, CV.GLM = cv, CV.MLM = cv,
               method = c("GLM", "MLM"), ncpus = 4, file.output = FALSE, verbose = FALSE)
    cat(sprintf("%s, covariates %s: %.0f s\n", tr, set, proc.time()[3] - t0))
    if (first) { cat("returned objects:", names(mvp), "\n"); first <- FALSE }
    for (meth in c("GLM", "MLM")) {
      res <- as.matrix(if (meth == "GLM") mvp$glm.results else mvp$mlm.results)
      label <- if (set == "none") c(GLM = "GLM_noPC", MLM = "MLM_K")[[meth]] else c(GLM = "GLM_PC3", MLM = "MLM_K_PC3")[[meth]]
      p <- res[, ncol(res)]
      tab <- data.frame(snp = map[, 1], chr = map[, 2], pos = map[, 3],
                        effect = signif(res[, 1], 5), se = signif(res[, 2], 5), p = signif(p, 5))
      write.table(tab, gzfile(paste0("results/gwas/", tr, "_", label, ".tsv.gz")),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      summ[[length(summ) + 1]] <- data.frame(trait = tr, model = label, lambda = round(lambda(p), 3),
          n_bonferroni = sum(p < bonf, na.rm = TRUE), n_1_over_m = sum(p < sugg, na.rm = TRUE),
          min_p = signif(min(p, na.rm = TRUE), 3))
      top <- head(tab[order(tab$p), ], 20)
      tops[[length(tops) + 1]] <- cbind(trait = tr, model = label, top)
    }
  }
}
summ <- do.call(rbind, summ); tops <- do.call(rbind, tops)
write.table(summ, "results/gwas_lambda.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tops, "results/gwas_top_hits.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n== lambda and significant SNP counts ==\n"); print(summ, row.names = FALSE)
cat("\n== top 5 SNPs per model ==\n")
for (k in split(tops, list(tops$trait, tops$model), drop = TRUE)) {
  cat("\n", k$trait[1], k$model[1], "\n"); print(head(k[, c("snp", "chr", "pos", "effect", "p")], 5), row.names = FALSE)
}
