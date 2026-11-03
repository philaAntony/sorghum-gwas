suppressMessages(library(data.table))
raw <- fread("data/work/pruned.raw")
idc <- intersect(names(raw), c("#FID","FID","IID","PAT","MAT","SEX","PHENOTYPE"))
M <- as.matrix(raw[, setdiff(names(raw), idc), with = FALSE]); rownames(M) <- raw$IID
sc <- read.delim("results/pca_scores.tsv")
sc <- sc[match(rownames(M), sc$taxa), ]
chr <- as.integer(sub(":.*", "", colnames(M)))
pos <- as.integer(sub("_.*", "", sub("^[0-9]+:", "", colnames(M))))
cat("example SNP names:", head(colnames(M), 2), "\n")
for (i in 1:6) {
  r <- abs(cor(M, sc[[paste0("PC", i)]]))
  top <- order(r, decreasing = TRUE)[1:round(0.01 * length(r))]
  tab <- sort(table(chr[top]), decreasing = TRUE)
  cat(sprintf("PC%d: top 1%% SNPs by chromosome: %s\n", i, paste(names(tab), tab, sep = "=", collapse = " ")))
  cat(sprintf("     strongest SNP: chr%d:%d (|r| = %.2f)\n", chr[top[1]], pos[top[1]], r[top[1]]))
}
