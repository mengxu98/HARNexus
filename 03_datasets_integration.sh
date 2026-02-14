#!/bin/bash

set -e

source "code/functions/utils.sh"

code_dir="code/datasets"
check_command Rscript

log_message "Start processing 27 datasets integration..."
if should_process "../../data/BrainData/integration/objects_list_raw.rds"; then
  run_r_script "$code_dir" "datasets_integration_01.R" "objects_list_raw"
else
  log_message "'objects_list_raw.rds' already processed!"
fi

if should_process "../../data/BrainData/integration/objects_list_processed.rds"; then
  run_r_script "$code_dir" "datasets_integration_02.R" "objects_list_processed"
else
  log_message "'objects_list_processed.rds' already processed!"
fi

if should_process "../../data/BrainData/integration/objects_raw.rds"; then
  run_r_script "$code_dir" "datasets_integration_02.R" "objects_raw"
else
  log_message "'objects_raw.rds' already processed!"
fi

if should_process "../../data/BrainData/integration/objects_filtered.rds"; then
  run_r_script "$code_dir" "datasets_integration_02.R" "objects_filtered"
else
  log_message "'objects_filtered.rds' already processed!"
fi

if should_process "../../data/BrainData/integration/objects_integrated.rds"; then
  run_r_script "$code_dir" "datasets_integration_03.R" "objects_integrated"
else
  log_message "'objects_integrated.rds' already processed!"
fi
