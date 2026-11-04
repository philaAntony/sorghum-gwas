#!/usr/bin/env bash
# 0/1/2 genotype export of the LD-pruned SNPs, read by the PCA and kinship script
set -euo pipefail
cd "$(dirname "$0")/.."
plink2 --bfile data/work/snps_pruned --chr-set 10 --export A --out data/work/pruned > /dev/null
ls -lh data/work/pruned.raw
