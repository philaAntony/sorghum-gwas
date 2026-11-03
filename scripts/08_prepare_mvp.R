suppressMessages(library(rMVP))
options(rMVP.OutputLog2File = FALSE)
# PLINK binary (MAF-filtered SNPs) -> rMVP format, with the VanRaden kinship matrix
MVP.Data(fileBed = "data/work/snps_maf", fileKin = TRUE, filePC = FALSE,
         out = "data/work/mvp_maf", ncpus = 2)
print(list.files("data/work", pattern = "^mvp_maf"))
