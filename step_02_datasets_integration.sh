#!/bin/bash

set -e

# Source the enhanced download utilities and log functions
source "code/functions/utils.sh"

code_dir="code/datasets"

log_message "Start processing 27 datasets integration..."
if should_process "../../data/BrainData/integration/objects_list_raw.rds"; then
  Rscript $code_dir/datasets_integration_01.R
  log_success "'objects_list_raw.rds' processed successfully!"
else
  log_message "'objects_list_raw.rds' already processed!"
fi

if should_process "../../data/BrainData/integration/objects_list_processed.rds"; then
  Rscript $code_dir/datasets_integration_02.R
  log_success "'objects_list_processed.rds' processed successfully!"
else
  log_message "'objects_list_processed.rds' already processed!"
fi

if should_process "../../data/BrainData/integration/objects_raw.rds"; then
  Rscript $code_dir/datasets_integration_02.R
  log_success "'objects_raw.rds' processed successfully!"
else
  log_message "'objects_raw.rds' already processed!"
fi

if should_process "../../data/BrainData/integration/objects_filtered.rds"; then
  Rscript $code_dir/datasets_integration_02.R
  log_success "'objects_filtered.rds' processed successfully!"
else
  log_message "'objects_filtered.rds' already processed!"
fi

if should_process "../../data/BrainData/integration/objects_integrated.rds"; then
  Rscript $code_dir/datasets_integration_03.R
  log_success "'objects_integrated.rds' processed successfully!"
else
  log_message "'objects_integrated.rds' already processed!"
fi
