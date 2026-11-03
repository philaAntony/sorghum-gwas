dir.create("results", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)
p <- read.delim("metadata/pheno_2021_HN.tsv", na.strings = "NA")
keep <- readLines("metadata/keep_samples.txt")
p <- p[p$taxa %in% keep, ]
cat("lines:", nrow(p), "\n")

skew <- function(x) { x <- x[!is.na(x)]; mean((x - mean(x))^3) / sd(x)^3 }
flag <- data.frame()
for (v in c("plant_height", "days_to_flower")) {
  x <- p[[v]]
  z <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  cat(sprintf("\n%s: n=%d mean=%.2f sd=%.2f median=%.2f min=%.2f max=%.2f skewness=%.2f\n",
      v, sum(!is.na(x)), mean(x, na.rm = TRUE), sd(x, na.rm = TRUE), median(x, na.rm = TRUE),
      min(x, na.rm = TRUE), max(x, na.rm = TRUE), skew(x)))
  if (v == "plant_height") cat(sprintf("skewness after log transform: %.2f\n", skew(log(x))))
  out <- which(abs(z) > 3)
  cat("lines with |z| > 3:", length(out), "\n")
  if (length(out)) flag <- rbind(flag, data.frame(trait = v, taxa = p$taxa[out], value = x[out], z = round(z[out], 2)))
}
write.table(flag, "results/phenotype_flagged.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(flag)
cat("\ncorrelation, height vs days to flower:", round(cor(p$plant_height, p$days_to_flower, use = "complete.obs"), 3), "\n")
cat("lines with a single height record:", sum(p$n_ph == 1, na.rm = TRUE), "\n")

png("figures/phenotype_distributions.png", width = 2000, height = 1600, res = 300)
par(mfrow = c(2, 2))
hist(p$plant_height, breaks = 40, main = "Plant height (m)", xlab = "line mean")
hist(log(p$plant_height), breaks = 40, main = "log plant height", xlab = "log of line mean")
hist(p$days_to_flower, breaks = 40, main = "Days to flower", xlab = "line mean")
qqnorm(log(p$plant_height), main = "QQ plot, log plant height"); qqline(log(p$plant_height))
dev.off()
