#!/usr/bin/env bash
# Runs the whole workflow from the repository root: bash scripts/run_all.sh
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data/work results figures metadata
run() { echo; echo "[$(date +%H:%M:%S)] $1"; shift; "$@"; }

run "00 check inputs"             bash scripts/00_check_inputs.sh
run "00 download annotation"      bash scripts/00_download_annotation.sh
run "01 phenotype table"          python3 scripts/01_make_phenotype_table.py
run "01b sample list"             bash scripts/01b_keep_samples.sh
run "02 genotype QC"              bash scripts/02_genotype_qc.sh
run "03 phenotype QC"             Rscript scripts/03_phenotype_qc.R
run "04 GWAS phenotype files"     Rscript scripts/04_make_gwas_phenotypes.R
run "05 line groups"              python3 scripts/05_line_groups.py
run "05b export pruned SNPs"      bash scripts/05b_export_pruned.sh
run "06 population structure"     Rscript scripts/06_population_structure.R
run "07 PC drivers"               Rscript scripts/07_pc_drivers.R
run "check rMVP"                  Rscript -e 'if (!requireNamespace("rMVP", quietly = TRUE)) stop("rMVP not installed: run bash scripts/install_rmvp.sh")'
run "08 prepare rMVP data"        Rscript scripts/08_prepare_mvp.R
run "09 GWAS"                     Rscript scripts/09_run_gwas.R
run "10 plots, loci, overlap"     Rscript scripts/10_plots.R
run "11 sensitivity analyses"     Rscript scripts/11_sensitivity.R
run "12 LD-based regions"         bash scripts/12_ld_regions.sh
run "12b lead SNP checks"         bash scripts/12b_lead_snp_checks.sh
run "13 chr6 carrier check"       Rscript scripts/13_chr6_check.R
run "14 annotation build check"   python3 scripts/14_build_check.py
run "15 candidate genes"          python3 scripts/15_candidate_genes.py
run "16 known genes"              python3 scripts/16_known_genes.py
echo; echo "[$(date +%H:%M:%S)] finished"
