p <- read.delim("metadata/pheno_2021_HN.tsv", na.strings = "NA")
keep <- readLines("metadata/keep_samples.txt")
p <- p[p$taxa %in% keep, ]
out <- data.frame(taxa = p$taxa,
                  log_plant_height = round(log(p$plant_height), 5),
                  days_to_flower = p$days_to_flower)
write.table(out, "metadata/pheno_gwas.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
# PLINK format (the .fam from script 02 uses family ID 0)
plink <- data.frame(`#FID` = 0, IID = out$taxa, log_plant_height = out$log_plant_height,
                    days_to_flower = out$days_to_flower, check.names = FALSE)
write.table(plink, "metadata/pheno_plink.txt", sep = "\t", quote = FALSE, row.names = FALSE)
cat("lines:", nrow(out), "| missing log height:", sum(is.na(out$log_plant_height)),
    "| missing days to flower:", sum(is.na(out$days_to_flower)), "\n")
