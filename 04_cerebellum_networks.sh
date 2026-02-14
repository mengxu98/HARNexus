#!/bin/bash

set -e

source "code/functions/utils.sh"

code_dir="code/hic"
res_dir="results/gse97942/"
if [ ! -d "$res_dir" ]; then
  mkdir -p "$res_dir"
fi

check_command Rscript

log_message "Start processing transcription factors..."
run_r_script "$code_dir" "00_tfs.R" "Transcription factors"

if [ ! -f "../../data/BrainData/processed/GSE97942/GSE97942_cerebellum_processed.rds" ]; then
  log_message "Start processing GSE97942 cerebellum data..."
  run_r_script "$code_dir" "01_GSE97942_annotation.R" "GSE97942 cerebellum"
else
  log_message "GSE97942 cerebellum data already processed!"
fi

export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export GOTO_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

if [ ! -f "$res_dir/astro_4methods_networks.rds" ]; then
  run_r_script "$code_dir" "02_network_GSE97942.R" "GSE97942 astrocyte networks"
else
  log_message "GSE97942 astrocyte networks already processed!"
fi
