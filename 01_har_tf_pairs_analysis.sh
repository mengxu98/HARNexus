#!/bin/bash

set -e

source "code/functions/utils.sh"

code_dir="code/har_tf_pairs_analysis"
check_command Rscript
check_command python3

log_message "Start processing HAR-TF pairs analysis..."
run_r_script "$code_dir" "01_genome_data.R" "Genome data"
run_python_script "$code_dir" "02_motif_scanning_human.py" "Motif scanning human"
run_python_script "$code_dir" "03_motif_scanning_chimp.py" "Motif scanning chimp"
run_r_script "$code_dir" "04_tfs_processing.R" "TFs processing"
run_python_script "$code_dir" "05_motif_scoring.py" "Motif scoring"
run_r_script "$code_dir" "06_motifs_processing.R" "Motifs processing"
