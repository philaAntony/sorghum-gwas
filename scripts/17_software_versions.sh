#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p results
{
  echo "python: $(python3 --version 2>&1)"
  echo "plink2: $(plink2 --version 2>&1 | head -n 1)"
  echo "bcftools: $(bcftools --version | head -n 1)"
  Rscript -e 'cat("R:", R.version.string, "\n"); for (p in c("rMVP", "bigmemory", "data.table")) cat(p, as.character(packageVersion(p)), "\n")'
} > results/software_versions.txt
cat results/software_versions.txt
