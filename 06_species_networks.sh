#!/bin/bash

set -e

source "code/functions/utils.sh"

log_message "Checking required commands..."
check_command Rscript
check_command python3
log_success "All required commands are installed"

if [ ! -f "../../data/BrainData/processed/GSE192774/GSE192774_processed.rds" ]; then
  log_message "Start processing GSE192774 data..."
  run_r_script "code/datasets" "GSE192774.R" "GSE192774 data"
else
  log_message "GSE192774 data already processed!"
fi

log_message "Stage 1: Network Construction and Comparison"
run_r_script "code/species_networks" "01_network_construction.R" "Network construction"
run_python_script "code/species_networks" "02_network_comparison.py" "Network comparison" "1"
run_python_script "code/species_networks" "03_genes_list.py" "Generate gene lists" "1"

log_message "Stage 2: ATAC Data Preprocessing"
run_r_script "code/species_networks" "04_atac_union_peaks.R" "ATAC union peaks"
run_r_script "code/species_networks" "05_atac_daccre.R" "DAcCRE analysis (Differential Accessibility CRE)"
run_r_script "code/species_networks" "06_atac_peak2gene.R" "Peak-to-gene association"

log_message "Stage 3: Evolution Analysis"
run_r_script "code/species_networks" "07_daccre_evolution.R" "DAcCRE-based evolution analysis"

log_success "All analysis completed!"
