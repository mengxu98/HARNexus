#!/bin/bash

# Enhanced download script for ROSMAP RDS data files
# Usage: ./download_rosmap_rds-data.sh

set -e

# Source the enhanced download utilities
source "$(dirname "$0")/../functions/utils.sh"

# https://personal.broadinstitute.org/cboix/ad427_data/Data/Processed_data/

# Create data directory if it doesn't exist
DATA_DIR="../../data/BrainData/raw/ROSMAP/processed_data"

log_message "Starting ROSMAP RDS data download..."

# Define download list with expected file sizes (approximate)
DOWNLOAD_LIST="
https://personal.broadinstitute.org/cboix/ad427_data/Data/Processed_data/Astrocytes.rds|Astrocytes.rds|0
https://personal.broadinstitute.org/cboix/ad427_data/Data/Processed_data/Excitatory_neurons_set1.rds|Excitatory_neurons_set1.rds|0
https://personal.broadinstitute.org/cboix/ad427_data/Data/Processed_data/Excitatory_neurons_set2.rds|Excitatory_neurons_set2.rds|0
https://personal.broadinstitute.org/cboix/ad427_data/Data/Processed_data/Excitatory_neurons_set3.rds|Excitatory_neurons_set3.rds|0
https://personal.broadinstitute.org/cboix/ad427_data/Data/Processed_data/Inhibitory_neurons.rds|Inhibitory_neurons.rds|0
https://personal.broadinstitute.org/cboix/ad427_data/Data/Processed_data/OPCs.rds|OPCs.rds|0
https://personal.broadinstitute.org/cboix/ad427_data/Data/Processed_data/Oligodendrocytes.rds|Oligodendrocytes.rds|0
https://personal.broadinstitute.org/cboix/ad427_data/Data/Processed_data/Vasculature_cells.rds|Vasculature_cells.rds|0
https://personal.broadinstitute.org/cboix/ad427_data/Data/Processed_data/Immune_cells.rds|Immune_cells.rds|0
"

# Perform batch download
batch_download "$DOWNLOAD_LIST" "$DATA_DIR" 5

# Clean up temporary files
cleanup_temp_files "$DATA_DIR"

log_success "ROSMAP RDS data download completed!"
