#!/bin/bash

# Download script for external data files
# Usage: ./download_data.sh

set -e

# https://cells.ucsc.edu/?ds=rosmap-ad-aging-brain+ad-atac+integration

# Create data directory if it doesn't exist
DATA_DIR="../../data/BrainData/ROSMAP/RNA+ATAC"
mkdir -p "$DATA_DIR"

echo "Downloading ROSMAP data..."

if [ ! -f "$DATA_DIR/matrix.mtx.gz" ]; then
  echo "Downloading matrix.mtx.gz..."
  curl -L -o "$DATA_DIR/matrix.mtx.gz" "https://cells.ucsc.edu/rosmap-ad-aging-brain/ad-atac/integration/matrix.mtx.gz"
fi

if [ ! -f "$DATA_DIR/features.tsv.gz" ]; then
  echo "Downloading features.tsv.gz..."
  curl -L -o "$DATA_DIR/features.tsv.gz" "https://cells.ucsc.edu/rosmap-ad-aging-brain/ad-atac/integration/features.tsv.gz"
fi

if [ ! -f "$DATA_DIR/barcodes.tsv.gz" ]; then
  echo "Downloading barcodes.tsv.gz..."
  curl -L -o "$DATA_DIR/barcodes.tsv.gz" "https://cells.ucsc.edu/rosmap-ad-aging-brain/ad-atac/integration/barcodes.tsv.gz"
fi

if [ ! -f "$DATA_DIR/meta.tsv" ]; then
  echo "Downloading meta.tsv..."
  curl -L -o "$DATA_DIR/meta.tsv" "https://cells.ucsc.edu/rosmap-ad-aging-brain/ad-atac/integration/meta.tsv"
fi

if [ ! -f "$DATA_DIR/UMAP_coordinates.coords.tsv.gz" ]; then
  echo "Downloading UMAP_coordinates.coords.tsv.gz..."
  curl -L -o "$DATA_DIR/UMAP_coordinates.coords.tsv.gz" "https://cells.ucsc.edu/rosmap-ad-aging-brain/ad-atac/integration/UMAP_coordinates.coords.tsv.gz"
fi

echo "Download completed successfully!"
