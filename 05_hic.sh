#!/bin/bash

set -e

source "code/functions/utils.sh"

code_dir="code/hic"
res_dir="results/hic/"
if [ ! -d "$res_dir" ]; then
  mkdir -p "$res_dir"
fi

check_command Rscript
check_command python3

if [ ! -f "results/hic/HAR_gene_HiC_supported.csv" ]; then
  log_message "Start processing Hi-C data..."
  run_python_script "$code_dir" "03_hic_data.py" "Hi-C data"
else
  log_message "Hi-C data already processed!"
fi

if [ ! -f "$res_dir/intersection_results.rds" ]; then
  run_r_script "$code_dir" "04_hic_intersect_analysis.R" "Hi-C intersect analysis"
else
  log_message "Hi-C intersect analysis already processed!"
fi
