#!/usr/bin/env bash
# Checks that the Dryad files needed by the pipeline are in data/raw/ (Dryad blocks command-line downloads).
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data/raw
VCF=data/raw/B.SDP_SAP_815_maf001_het10_imputed.recode.vcf.gz
ZIP=data/raw/2.SDP_SAP_f2021_field_phenotype.zip
if [ ! -f "$VCF" ] || { [ ! -f "$ZIP" ] && ! ls data/raw/2.*/*.csv >/dev/null 2>&1; }; then
  cat << 'MSG'
Missing input data. Download these two files in a browser from the Dryad dataset
https://doi.org/10.5061/dryad.3j9kd51w9 and place them in data/raw/:
  B.SDP_SAP_815_maf001_het10_imputed.recode.vcf.gz
  2.SDP_SAP_f2021_field_phenotype.zip
(Dryad refuses command-line downloads.)
MSG
  exit 1
fi
if ! ls data/raw/2.*/*.csv >/dev/null 2>&1; then
  python3 -c "import zipfile; zipfile.ZipFile('$ZIP').extractall('data/raw/2.SDP_SAP_f2021_field_phenotype')"
fi
gzip -t "$VCF" && echo "inputs OK"
