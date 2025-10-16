#!/bin/bash

# Enhanced download script for ROSMAP UCSC snRNA-seq data
# Usage: ./download_rosmap_ucsc_snRNA-seq.sh

set -e

# Source the enhanced download utilities
source "$(dirname "$0")/download_utils.sh"

# https://cells.ucsc.edu/?ds=rosmap-ad-aging-brain+ad-atac+rna-hq

# Create data directory if it doesn't exist
DATA_DIR="../../data/BrainData/ROSMAP/RNA"

log_message "Starting ROSMAP UCSC snRNA-seq data download..."

# Define download list
DOWNLOAD_LIST="
https://cells.ucsc.edu/rosmap-ad-aging-brain/ad-atac/rna-hq/matrix.mtx.gz|matrix.mtx.gz|0
https://cells.ucsc.edu/rosmap-ad-aging-brain/ad-atac/rna-hq/features.tsv.gz|features.tsv.gz|0
https://cells.ucsc.edu/rosmap-ad-aging-brain/ad-atac/rna-hq/barcodes.tsv.gz|barcodes.tsv.gz|0
https://cells.ucsc.edu/rosmap-ad-aging-brain/ad-atac/rna-hq/meta.tsv|meta.tsv|0
https://cells.ucsc.edu/rosmap-ad-aging-brain/ad-atac/rna-hq/UMAP_coordinates.coords.tsv.gz|UMAP_coordinates.coords.tsv.gz|0
"

# Perform batch download
batch_download "$DOWNLOAD_LIST" "$DATA_DIR" 5

# Clean up temporary files
cleanup_temp_files "$DATA_DIR"

log_success "ROSMAP UCSC snRNA-seq data download completed!"
