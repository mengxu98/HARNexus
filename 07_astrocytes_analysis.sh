#!/bin/bash

set -e

source "code/functions/utils.sh"

log_message "Checking required commands..."
check_command Rscript
check_command python3
log_message "All required commands are installed" --message-type success

log_message "Start processing astrocytes analysis..."
run_python_script "code/networks" "08_pfc_celltype_targets_analysis.py" "PFC cell type targets analysis"
run_r_script "code/networks" "09_region_object.R" "Region object"
run_r_script "code/networks" "10_pfc_astrocytes_analysis.R" "PFC astrocytes analysis"
