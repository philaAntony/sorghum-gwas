K <- readRDS("data/work/kinship_vanraden.rds")
R <- K / sqrt(outer(diag(K), diag(K)))
d <- read.table("data/work/chr6_alt_dosage.txt", col.names = c("taxa", "dose"), colClasses = c("character", "numeric"))
g <- read.delim("metadata/line_groups.tsv", colClasses = "character")
ph <- read.delim("metadata/pheno_gwas.tsv", colClasses = c("character", "numeric", "numeric"))
car <- d$taxa[d$dose > 0]
car <- car[car %in% rownames(R)]
mean_rel <- function(ids) { m <- R[ids, ids]; mean(m[upper.tri(m)]) }
set.seed(1)
sap <- intersect(g$taxa[g$panel == "SAP"], rownames(R))
r_all <- replicate(1000, mean_rel(sample(rownames(R), length(car))))
r_sap <- replicate(1000, mean_rel(sample(sap, length(car))))
cat(sprintf("carriers: %d | mean pairwise relatedness among carriers: %.3f\n", length(car), mean_rel(car)))
cat(sprintf("random %d lines: mean %.3f (95%% range %.3f to %.3f)\n", length(car), mean(r_all), quantile(r_all, .025), quantile(r_all, .975)))
cat(sprintf("random %d SAP lines: mean %.3f (95%% range %.3f to %.3f)\n", length(car), mean(r_sap), quantile(r_sap, .025), quantile(r_sap, .975)))
m <- merge(ph, g, by = "taxa"); m$carrier <- m$taxa %in% car
cat("\ndays to flower by panel and carrier status:\n")
print(aggregate(days_to_flower ~ panel + carrier, m, function(x) c(n = length(x), mean = round(mean(x), 2))))
s <- m[m$panel == "SAP", ]
cat("\nwithin SAP, carrier vs non-carrier t-test p =", signif(t.test(days_to_flower ~ carrier, s)$p.value, 3), "\n")
