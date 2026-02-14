#!/bin/bash

set -e

source "code/functions/utils.sh"

log_message "Building Supplementary Tables (xlsx and CSVs)"
echo ""

code_dir="code/supplementary_file"
res_dir="results/supplementary_file"
if [ ! -d "$res_dir" ]; then
  mkdir -p "$res_dir"
fi

check_command python3

run_python_script "$code_dir" "prepare_supplementary_tables.py" "Supplementary tables"

log_success "Supplementary tables completed!"
echo ""
