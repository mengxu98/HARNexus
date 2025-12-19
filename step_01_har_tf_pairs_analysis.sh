#!/bin/bash

set -e

# Source the enhanced download utilities and log functions
source "code/functions/utils.sh"

code_dir="code/har_tf_pairs_analysis"

log_message "Start processing HAR-TF pairs analysis..."
Rscript $code_dir/step01a_genome_data.R
python $code_dir/step01b_motif_scanning_chimp
python $code_dir/step01b_motif_scanning_human.py
Rscript $code_dir/step01c_tfs_processing
python $code_dir/step01d_motif_scoring.py
Rscript $code_dir/step01e_motif_aggregation.R
