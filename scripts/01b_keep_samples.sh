#!/usr/bin/env bash
# Lines with both traits, excluding the non-accession "Border"
set -euo pipefail
cd "$(dirname "$0")/.."
awk -F'\t' 'NR>1 && $2!="NA" && $4!="NA" && $1!="Border" {print $1}' metadata/pheno_2021_HN.tsv > metadata/keep_samples.txt
echo "lines kept: $(wc -l < metadata/keep_samples.txt)"
