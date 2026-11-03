suppressMessages(library(data.table))
dir.create("results", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)

# ---- genotypes (0/1/2, 23,900 pruned SNPs) ----
raw <- fread("data/work/pruned.raw")
idcols <- intersect(names(raw), c("#FID", "FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"))
M <- as.matrix(raw[, setdiff(names(raw), idcols), with = FALSE])
rownames(M) <- raw$IID
cat("genotype matrix:", dim(M), "| missing values:", sum(is.na(M)), "\n")
stopifnot(!anyNA(M))

# ---- VanRaden kinship: centre by 2p, scale by 2*sum(p(1-p)) ----
p <- colMeans(M) / 2
W <- sweep(M, 2, 2 * p)
K <- tcrossprod(W) / (2 * sum(p * (1 - p)))
d <- diag(K)
off <- K[upper.tri(K)]
cat("\nKinship diagonal: mean", round(mean(d), 3), "range", round(range(d), 3), "\n")
cat("Off-diagonal quantiles:\n"); print(round(quantile(off, c(0.01, 0.25, 0.5, 0.75, 0.99)), 3))

# standardised relatedness between lines: K_ij / sqrt(K_ii K_jj)
R <- K / sqrt(outer(d, d))
idx <- which(R > 0.9 & upper.tri(R), arr.ind = TRUE)
dups <- data.frame(line1 = rownames(R)[idx[, 1]], line2 = colnames(R)[idx[, 2]], relatedness = round(R[idx], 3))
dups <- dups[order(-dups$relatedness), ]
write.table(dups, "results/near_duplicates.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nPairs with standardised relatedness > 0.9:", nrow(dups), "\n"); print(head(dups, 10))

# ---- PCA from the eigen-decomposition of K ----
e <- eigen(K, symmetric = TRUE)
pct <- 100 * e$values / sum(e$values)
cat("\nVariance explained by PC1-10 (%):", round(pct[1:10], 2), "\n")
cat("Cumulative PC1-5:", round(sum(pct[1:5]), 1), "%\n")
pcs <- sweep(e$vectors[, 1:10], 2, sqrt(pmax(e$values[1:10], 0)), "*")
colnames(pcs) <- paste0("PC", 1:10); rownames(pcs) <- rownames(K)
write.table(data.frame(PC = 1:20, percent = round(pct[1:20], 3)), "results/pca_variance.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- merge line information ----
g <- read.delim("metadata/line_groups.tsv", stringsAsFactors = FALSE)
ph <- read.delim("metadata/pheno_gwas.tsv", stringsAsFactors = FALSE)
df <- data.frame(taxa = rownames(pcs), pcs, stringsAsFactors = FALSE)
df <- merge(merge(df, g, by = "taxa", all.x = TRUE), ph, by = "taxa", all.x = TRUE)
set.seed(1)
km <- kmeans(pcs[, 1:5], centers = 4, nstart = 25)
df$cluster <- km$cluster[match(df$taxa, rownames(pcs))]
write.table(df, "results/pca_scores.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nk-means clusters (k = 4, exploratory):\n"); print(table(df$cluster))
cat("\ncluster x panel:\n"); print(table(df$cluster, df$panel))
cat("\ncluster x name type:\n"); print(table(df$cluster, df$name_type))

# ---- structure versus traits ----
cat("\nTraits versus structure:\n")
for (v in c("log_plant_height", "days_to_flower")) {
  fit <- lm(df[[v]] ~ PC1 + PC2 + PC3 + PC4 + PC5, data = df)
  cat(sprintf("  %s: R2 explained by PC1-5 = %.3f; correlations with PC1..5: %s\n", v,
      summary(fit)$r.squared,
      paste(round(sapply(1:5, function(i) cor(df[[v]], df[[paste0("PC", i)]], use = "complete.obs")), 2), collapse = " ")))
  cat("  means by panel:", paste(names(tapply(df[[v]], df$panel, mean, na.rm = TRUE)),
      round(tapply(df[[v]], df$panel, mean, na.rm = TRUE), 3), collapse = "; "), "\n")
  cat("  means by cluster:", paste(round(tapply(df[[v]], df$cluster, mean, na.rm = TRUE), 3), collapse = " "), "\n")
}

# ---- the most heterozygous lines ----
het <- read.delim("results/sample_het.tsv", header = FALSE, col.names = c("taxa", "het"))
top <- head(het$taxa[order(-het$het)], 5)
Rn <- R; diag(Rn) <- NA
cat("\nHighest-heterozygosity lines:\n")
for (t in top) {
  r <- df[df$taxa == t, ]
  cat(sprintf("  %s: het %.3f | PC1 %.1f PC2 %.1f | panel %s | cluster %d | max relatedness to another line %.2f\n",
      t, het$het[het$taxa == t], r$PC1, r$PC2, r$panel, r$cluster, max(Rn[t, ], na.rm = TRUE)))
}

# ---- figures ----
lab <- function(i) sprintf("PC%d (%.1f%%)", i, pct[i])
pal_panel <- c("SAP" = "#1b9e77", "non-SAP" = "#d95f02")
pal_type <- c(PI = "grey55", SC = "#1f78b4", BTx_RTx = "#e31a1c", other = "#33a02c")
pal_clu <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")

png("figures/pca_groups.png", width = 3300, height = 1200, res = 300)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
plot(df$PC1, df$PC2, col = pal_panel[df$panel], pch = 19, cex = 0.6, xlab = lab(1), ylab = lab(2), main = "SAP membership")
legend("topright", names(pal_panel), col = pal_panel, pch = 19, bty = "n", cex = 0.8)
plot(df$PC1, df$PC2, col = pal_type[df$name_type], pch = 19, cex = 0.6, xlab = lab(1), ylab = lab(2), main = "Line name type")
legend("topright", names(pal_type), col = pal_type, pch = 19, bty = "n", cex = 0.8)
plot(df$PC1, df$PC2, col = pal_clu[df$cluster], pch = 19, cex = 0.6, xlab = lab(1), ylab = lab(2), main = "k-means clusters (k=4)")
dev.off()

ramp <- colorRampPalette(c("#2166ac", "grey85", "#b2182b"))(50)
png("figures/pca_traits.png", width = 2400, height = 1200, res = 300)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
for (v in c("log_plant_height", "days_to_flower")) {
  cc <- ramp[cut(df[[v]], 50, labels = FALSE)]
  plot(df$PC1, df$PC2, col = cc, pch = 19, cex = 0.6, xlab = lab(1), ylab = lab(2), main = v)
}
dev.off()

png("figures/pca_scree.png", width = 1200, height = 1000, res = 300)
barplot(pct[1:20], names.arg = 1:20, xlab = "PC", ylab = "% variance", main = "Scree plot")
dev.off()

hc <- hclust(as.dist(1 - R), method = "average"); ord <- hc$order
n <- nrow(R)
png("figures/kinship_heatmap.png", width = 1600, height = 1600, res = 200)
image(1:n, 1:n, pmax(R[ord, ord], 0), col = colorRampPalette(c("white", "#fdae61", "#a50026"))(50),
      axes = FALSE, xlab = "lines (clustered)", ylab = "lines (clustered)", main = "Standardised relatedness")
dev.off()
saveRDS(K, "data/work/kinship_vanraden.rds")
cat("\nfigures written\n")
