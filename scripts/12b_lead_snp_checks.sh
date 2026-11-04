#!/usr/bin/env bash
# LD between the two chr9 lead SNPs, and the ALT-allele dosage of the chr6 flowering lead SNP
set -euo pipefail
cd "$(dirname "$0")/.."
ht9=$(awk -F'\t' '$1=="log_plant_height" && $2==9 {print $6}' results/gwas_loci.tsv)
fl9=$(awk -F'\t' '$1=="days_to_flower" && $2==9 {print $6}' results/gwas_loci.tsv)
fl6=$(awk -F'\t' '$1=="days_to_flower" && $2==6 {print $6}' results/gwas_loci.tsv)
echo "lead SNPs: height chr9 $ht9 | flowering chr9 $fl9 | flowering chr6 $fl6"
plink2 --bfile data/work/snps_maf --chr-set 10 --r2-unphased --ld-snp "$ht9" \
       --ld-window-kb 3000 --ld-window-r2 0 --out data/work/ld_lead > /dev/null
{ head -n 1 data/work/ld_lead.vcor; grep -F "$fl9" data/work/ld_lead.vcor; } > results/ld_chr9_leads.tsv
plink2 --bfile data/work/snps_maf --chr-set 10 --snp "$fl6" --export A --out data/work/lead_chr6 > /dev/null
alt=$(awk -v s="$fl6" '$2==s {print $4}' data/work/qc.afreq)
counted=$(head -n 1 data/work/lead_chr6.raw | tr '\t' '\n' | tail -n 1 | sed 's/.*_//')
awk -v alt="$alt" -v cnt="$counted" 'NR>1 { d = (cnt==alt) ? $7 : 2-$7; print $2, d }' data/work/lead_chr6.raw > data/work/chr6_alt_dosage.txt
echo "ALT $alt, counted $counted; dosage table: $(wc -l < data/work/chr6_alt_dosage.txt) lines"
