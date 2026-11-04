#!/usr/bin/env bash
# rMVP is a CRAN package (needs a C++ compiler); installed into the conda environment
Rscript -e 'install.packages("rMVP", repos = c("https://cran.r-universe.dev", "https://cloud.r-project.org"))'
