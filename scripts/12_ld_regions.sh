#!/usr/bin/env bash
# Candidate regions: significant SNPs merged with SNPs in LD (r2 >= 0.5) with each lead SNP
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p results/candidates
printf "trait\tchr\tlead_snp\tsig_start\tsig_end\tn_ld_snps\tregion_start\tregion_end\tregion_mb\n" > results/candidate_regions.tsv
tail -n +2 results/gwas_loci.tsv | while IFS=$'\t' read -r trait chr start end n lead leadpos leadp eff; do
  out=results/candidates/ld_${trait}_chr${chr}
  plink2 --bfile data/work/snps_maf --chr-set 10 --r2-unphased --ld-snp "$lead" \
         --ld-window-kb 3000 --ld-window-r2 0.5 --out "$out" > /dev/null
  awk -F'\t' -v t="$trait" -v c="$chr" -v l="$lead" -v s="$start" -v e="$end" '
    BEGIN { lo = s; hi = e; n = 0 }
    FNR > 1 { n++; if ($5 < lo) lo = $5; if ($5 > hi) hi = $5 }
    END { printf "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%.2f\n", t, c, l, s, e, n, lo, hi, (hi - lo) / 1e6 }' \
    "$out.vcor" >> results/candidate_regions.tsv
done
column -t results/candidate_regions.tsv
