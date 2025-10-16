#!/bin/bash

# Download script for external data files
# Usage: ./download_data.sh

set -e

# https://personal.broadinstitute.org/cboix/ad427_data/Data/Processed_data/

# Create data directory if it doesn't exist
DATA_DIR="../../data/BrainData/ROSMAP/processed_data"
mkdir -p "$DATA_DIR"


echo "Downloading ROSMAP data..."

if [ ! -f "$DATA_DIR/Astrocytes.rds" ]; then
  echo "Downloading Astrocytes.rds..."
  curl -L -o "$DATA_DIR/Astrocytes.rds" "https://personal.broadinstitute.org/cboix/ad427_data/Data/Processed_data/Astrocytes.rds"
fi

if [ ! -f "$DATA_DIR/Excitatory_neurons_set1.rds" ]; then
  echo "Downloading Excitatory_neurons_set1.rds..."
  curl -L -o "$DATA_DIR/Excitatory_neurons_set1.rds" "https://personal.broadinstitute.org/cboix/ad427_data/Data/Processed_data/Excitatory_neurons_set1.rds"
fi

if [ ! -f "$DATA_DIR/Excitatory_neurons_set2.rds" ]; then
  echo "Downloading Excitatory_neurons_set2.rds..."
  curl -L -o "$DATA_DIR/Excitatory_neurons_set2.rds" "https://personal.broadinstitute.org/cboix/ad427_data/Data/Processed_data/Excitatory_neurons_set2.rds"
fi

if [ ! -f "$DATA_DIR/Excitatory_neurons_set3.rds" ]; then
  echo "Downloading Excitatory_neurons_set3.rds..."
  curl -L -o "$DATA_DIR/Excitatory_neurons_set3.rds" "https://personal.broadinstitute.org/cboix/ad427_data/Data/Processed_data/Excitatory_neurons_set3.rds"
fi

if [ ! -f "$DATA_DIR/Inhibitory_neurons.rds" ]; then
  echo "Downloading Inhibitory_neurons.rds..."
  curl -L -o "$DATA_DIR/Inhibitory_neurons.rds" "https://personal.broadinstitute.org/cboix/ad427_data/Data/Processed_data/Inhibitory_neurons.rds"
fi

if [ ! -f "$DATA_DIR/OPCs.rds" ]; then
  echo "Downloading OPCs.rds..."
  curl -L -o "$DATA_DIR/OPCs.rds" "https://personal.broadinstitute.org/cboix/ad427_data/Data/Processed_data/OPCs.rds"
fi

if [ ! -f "$DATA_DIR/Oligodendrocytes.rds" ]; then
  echo "Downloading Oligodendrocytes.rds..."
  curl -L -o "$DATA_DIR/Oligodendrocytes.rds" "https://personal.broadinstitute.org/cboix/ad427_data/Data/Processed_data/Oligodendrocytes.rds"
fi

if [ ! -f "$DATA_DIR/Vasculature_cells.rds" ]; then
  echo "Downloading Vasculature_cells.rds..."
  curl -L -o "$DATA_DIR/Vasculature_cells.rds" "https://personal.broadinstitute.org/cboix/ad427_data/Data/Processed_data/Vasculature_cells.rds"
fi

if [ ! -f "$DATA_DIR/Immune_cells.rds" ]; then
  echo "Downloading Immune_cells.rds..."
  curl -L -o "$DATA_DIR/Immune_cells.rds" "https://personal.broadinstitute.org/cboix/ad427_data/Data/Processed_data/Immune_cells.rds"
fi

echo "Download completed successfully!"
