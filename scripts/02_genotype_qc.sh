#!/usr/bin/env bash
# Genotype QC. Run: bash scripts/02_genotype_qc.sh   (optional: MAF=0.01 bash ...)
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data/work results

V=data/raw/B.SDP_SAP_815_maf001_het10_imputed.recode.vcf.gz
MAF=${MAF:-0.05}

for i in 01 02 03 04 05 06 07 08 09 10; do echo "Chr$i ${i#0}"; done > metadata/chr_map.txt
n_raw_vars=$(bcftools view -H "$V" | wc -l)
n_raw_samp=$(bcftools query -l "$V" | wc -l)

# 1. biallelic SNPs in the kept samples (PLINK wants chromosomes 1-10)
bcftools view -S metadata/keep_samples.txt -m2 -M2 -v snps "$V" \
  | bcftools annotate --rename-chrs metadata/chr_map.txt -Oz -o data/work/snps_kept.vcf.gz
bcftools index -f data/work/snps_kept.vcf.gz
plink2 --vcf data/work/snps_kept.vcf.gz --chr-set 10 --set-all-var-ids '@:#' --make-bed --out data/work/snps_kept
n_snps=$(wc -l < data/work/snps_kept.bim)
n_samp=$(wc -l < data/work/snps_kept.fam)

# 2. statistics
plink2 --bfile data/work/snps_kept --chr-set 10 --freq --missing --sample-counts --out data/work/qc

# 3. MAF filter
plink2 --bfile data/work/snps_kept --chr-set 10 --maf "$MAF" --make-bed --out data/work/snps_maf
n_maf=$(wc -l < data/work/snps_maf.bim)

# 4. LD thinning (for PCA and kinship)
plink2 --bfile data/work/snps_maf --chr-set 10 --indep-pairwise 50 5 0.2 --out data/work/ld
n_ld=$(wc -l < data/work/ld.prune.in)
plink2 --bfile data/work/snps_maf --chr-set 10 --extract data/work/ld.prune.in --make-bed --out data/work/snps_pruned

# 5. tables
{
printf "stage\tsamples\tvariants\n"
printf "raw VCF\t%s\t%s\n" "$n_raw_samp" "$n_raw_vars"
printf "biallelic SNPs, kept samples\t%s\t%s\n" "$n_samp" "$n_snps"
printf "after MAF >= %s\t%s\t%s\n" "$MAF" "$n_samp" "$n_maf"
printf "after LD thinning (r2 < 0.2, 50-SNP window)\t%s\t%s\n" "$n_samp" "$n_ld"
} > results/qc_snp_counts.tsv

awk -F'\t' 'NR>1 {n=$3+$4+$5; printf "%s\t%.5f\n", $2, $5/n}' data/work/qc.scount > results/sample_het.tsv
{
  column -t -s$'\t' results/qc_snp_counts.tsv
  echo
  awk -F'\t' 'NR>1{ts+=$6; tv+=$7} END{printf "Ts/Tv ratio: %.3f\n", ts/tv}' data/work/qc.scount
  awk -F'\t' 'NR>1 && $12>0 {m++} END{print "samples with missing genotypes:", m+0}' data/work/qc.scount
  awk 'NR>1 {f=$6; m=(f<1-f)?f:1-f; if(m<0.01)a++; else if(m<0.05)b++; else if(m<0.1)c++; else d++}
       END{print "MAF bins  <0.01:", a+0, " 0.01-0.05:", b+0, " 0.05-0.1:", c+0, " >=0.1:", d+0}' data/work/qc.afreq
  echo "per-sample heterozygosity (min, median, max):"
  sort -t$'\t' -k2,2g results/sample_het.tsv | awk -F'\t' '{a[NR]=$2} END{print a[1], a[int(NR/2)], a[NR]}'
  echo "highest-heterozygosity samples:"
  sort -t$'\t' -k2,2gr results/sample_het.tsv | head -n 5
} | tee results/qc_summary.txt
