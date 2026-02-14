#!/bin/bash

set -e

source "code/functions/utils.sh"

code_dir="code/networks"
res_dir="results/networks/gse97942/"
if [ ! -d "$res_dir" ]; then
  mkdir -p "$res_dir"
fi

check_command Rscript

log_message "Start processing transcription factors..."
run_r_script "code/hic" "00_tfs.R" "Transcription factors"

log_message "Start processing all datasets data..."
run_r_script "$code_dir" "01_har_csn_data.R" "All datasets data"

# setting up the environment for sugon server
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export GOTO_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

log_message "Start processing HAR-CSN atlas networks..."
run_r_script "$code_dir" "02_har_csn_atlas.R" "HAR-CSN atlas networks"
