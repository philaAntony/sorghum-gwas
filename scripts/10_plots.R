suppressMessages(library(data.table))
dir.create("figures", showWarnings = FALSE)
traits <- c(log_plant_height = "log plant height", days_to_flower = "days to flower")
models <- c(GLM_noPC = "GLM, no PCs", GLM_PC3 = "GLM + PC1-3", MLM_K = "MLM, kinship", MLM_K_PC3 = "MLM, kinship + PC1-3")
lam <- read.delim("results/gwas_lambda.tsv")
lam_of <- function(tr, md) lam$lambda[lam$trait == tr & lam$model == md]

tabs <- list()
for (tr in names(traits)) for (md in names(models))
  tabs[[paste(tr, md)]] <- fread(cmd = paste("zcat", sprintf("results/gwas/%s_%s.tsv.gz", tr, md)))
m <- nrow(tabs[[1]])
bonf <- -log10(0.05 / m); sugg <- -log10(1 / m)
n_pruned <- length(readLines("data/work/ld.prune.in"))
cat(sprintf("m = %d; Bonferroni -log10 = %.2f; 1/m -log10 = %.2f; pruned-count Bonferroni -log10 = %.2f\n",
            m, bonf, sugg, -log10(0.05 / n_pruned)))

ref <- tabs[[1]]
chr_len <- tapply(ref$pos, ref$chr, max)
off <- c(0, cumsum(as.numeric(chr_len)))[seq_along(chr_len)]; names(off) <- names(chr_len)
mid <- off + chr_len / 2

manhattan <- function(d, ylim, main) {
  x <- d$pos + off[as.character(d$chr)]
  y <- -log10(pmax(d$p, 1e-300))
  col <- ifelse(d$chr %% 2 == 1, "grey35", "#2c7fb8")
  plot(x, y, pch = 19, cex = 0.35, col = col, xaxt = "n", xlab = "chromosome",
       ylab = expression(-log[10](p)), ylim = ylim, main = main, cex.main = 0.95)
  axis(1, at = mid, labels = names(mid), tick = FALSE, cex.axis = 0.8)
  abline(h = bonf, col = "red", lty = 2); abline(h = sugg, col = "blue", lty = 3)
}

for (tr in names(traits)) {
  ymax <- max(sapply(names(models), function(md) max(-log10(pmax(tabs[[paste(tr, md)]]$p, 1e-300)), na.rm = TRUE)))
  png(sprintf("figures/manhattan_%s.png", tr), width = 3200, height = 2600, res = 300)
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  for (md in names(models))
    manhattan(tabs[[paste(tr, md)]], c(0, ymax * 1.05),
              sprintf("%s: %s (lambda = %.2f)", traits[[tr]], models[[md]], lam_of(tr, md)))
  dev.off()
}

cols <- c(GLM_noPC = "#d95f02", GLM_PC3 = "#e7298a", MLM_K = "#1b9e77", MLM_K_PC3 = "#1f78b4")
png("figures/qq_plots.png", width = 2800, height = 1400, res = 300)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
for (tr in names(traits)) {
  n <- m; e <- -log10((seq_len(n) - 0.5) / n)
  allmax <- max(sapply(names(models), function(md) max(-log10(pmax(tabs[[paste(tr, md)]]$p, 1e-300)), na.rm = TRUE)))
  plot(NA, xlim = c(0, max(e)), ylim = c(0, allmax * 1.05), xlab = expression(Expected~-log[10](p)),
       ylab = expression(Observed~-log[10](p)), main = traits[[tr]])
  abline(0, 1, col = "black")
  for (md in names(models)) {
    o <- -log10(pmax(sort(tabs[[paste(tr, md)]]$p), 1e-300))
    points(e, o, pch = 19, cex = 0.3, col = cols[[md]])
  }
  legend("topleft", sprintf("%s (lambda %.2f)", models, sapply(names(models), function(md) lam_of(tr, md))),
         col = cols, pch = 19, bty = "n", cex = 0.7)
}
dev.off()

# loci: significant SNPs in the primary model, merged when within 1 Mb
loci <- function(d, thr = 0.05 / m, gap = 1e6) {
  s <- d[!is.na(p) & p < thr][order(chr, pos)]
  if (!nrow(s)) return(NULL)
  s[, locus := cumsum(c(TRUE, diff(pos) > gap | diff(chr) != 0))]
  s[, .(chr = chr[1], start = min(pos), end = max(pos), n_snps = .N,
        lead_snp = snp[which.min(p)], lead_pos = pos[which.min(p)],
        lead_p = min(p), lead_effect = effect[which.min(p)]), by = locus]
}
out <- list()
for (tr in names(traits)) { l <- loci(tabs[[paste(tr, "MLM_K_PC3")]]); if (!is.null(l)) out[[tr]] <- cbind(trait = tr, l) }
lc <- rbindlist(out); lc[, locus := NULL]
write.table(lc, "results/gwas_loci.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n== loci (MLM kinship + PC1-3, Bonferroni, SNPs merged within 1 Mb) ==\n"); print(lc)

# overlap of significant SNPs between models
ov <- list()
for (tr in names(traits)) {
  prim <- tabs[[paste(tr, "MLM_K_PC3")]]; ps <- prim$snp[!is.na(prim$p) & prim$p < 0.05 / m]
  for (md in names(models)) {
    d <- tabs[[paste(tr, md)]]; s <- d$snp[!is.na(d$p) & d$p < 0.05 / m]
    ov[[length(ov) + 1]] <- data.frame(trait = tr, model = md, n_sig = length(s),
      shared_with_primary = length(intersect(s, ps)), only_in_this_model = length(setdiff(s, ps)))
  }
}
ov <- do.call(rbind, ov)
write.table(ov, "results/gwas_model_overlap.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n== overlap of Bonferroni-significant SNPs with the primary model ==\n"); print(ov, row.names = FALSE)
