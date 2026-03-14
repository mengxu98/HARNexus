#!/bin/bash

set -e

source "code/functions/utils.sh"

code_dir="code/networks"

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


log_message "Start processing atlas processing..."
run_python_script "$code_dir" "03_atlas_processing.py" "Atlas processing"
run_python_script "$code_dir" "04_altas_statistics.py" "Statistics analysis"
run_python_script "$code_dir" "05_similarity_analysis.py" "Similarity analysis"
run_python_script "$code_dir" "06_all_stages_tf_target_pairs.py" "All stages TF-target pairs"
run_python_script "$code_dir" "07_all_stages_target_genes.py" "All stages target genes"
