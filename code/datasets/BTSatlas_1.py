#!/usr/bin/env python3
"""
Convert BTS_atlas.h5ad to a format compatible with scop's adata_to_srt function
and save count matrix and metadata as CSV files
paper: https://doi.org/10.1038/s12276-024-01328-6
data: https://zenodo.org/records/10939707
"""

import scanpy as sc
import pandas as pd  # Needed for DataFrame operations and to_csv
from scipy import sparse
import os
import sys

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from functions.utils import log_message

# Set paths
data_dir = "../../data/BrainData/raw/BTSatlas"
res_dir = "../../data/BrainData/processed/BTSatlas/"
os.makedirs(res_dir, exist_ok=True)

log_message("Start loading data...", message_type="running")

# Load h5ad file
h5ad_path = os.path.join(data_dir, "BTS_atlas.h5ad")
if not os.path.exists(h5ad_path):
    log_message("Error: {.path ", h5ad_path, "} does not exist", message_type="error")
    sys.exit(1)

log_message("Loading h5ad file: {.path ", h5ad_path, "}")
# Load without backed mode to ensure all data is accessible
# For large files, this will use more memory but ensures compatibility
adata = sc.read_h5ad(h5ad_path)

log_message("AnnData shape: {.val ", adata.shape, "}")
log_message("Number of cells: {.val ", adata.n_obs, "}")
log_message("Number of genes: {.val ", adata.n_vars, "}")

# Check available layers
if hasattr(adata, "layers") and len(adata.layers) > 0:
    log_message(
        "Available layers: {.val ", list(adata.layers.keys()), "}", message_type="info"
    )

# Extract count matrix - try multiple strategies
count_matrix = None
gene_names = None
var_df = None

# Strategy 1: Try raw layer first
if adata.raw is not None:
    try:
        log_message("Trying raw layer for count matrix...", message_type="info")
        count_matrix = adata.raw.X
        gene_names = adata.raw.var_names
        var_df = adata.raw.var.copy()
        log_message("Using raw layer for count matrix", message_type="success")
    except Exception as e:
        log_message(
            "Failed to access raw layer: {.val ", str(e), "}", message_type="warning"
        )

# Strategy 2: Try main X matrix
if count_matrix is None:
    try:
        log_message("Trying main X matrix...", message_type="info")
        count_matrix = adata.X
        gene_names = adata.var_names
        var_df = adata.var.copy()
        log_message("Using main X matrix for count matrix", message_type="success")
    except Exception as e:
        log_message(
            "Failed to access main X matrix: {.val ",
            str(e),
            "}",
            message_type="warning",
        )

# Strategy 3: Try layers if X is not available
if count_matrix is None and hasattr(adata, "layers") and len(adata.layers) > 0:
    # Try common layer names for count data
    layer_candidates = ["counts", "raw_counts", "logcounts"]
    for layer_name in layer_candidates:
        if layer_name in adata.layers:
            try:
                log_message(
                    "Trying layer '{.val ",
                    layer_name,
                    "}' for count matrix...",
                    message_type="info",
                )
                count_matrix = adata.layers[layer_name]
                gene_names = adata.var_names
                var_df = adata.var.copy()
                log_message(
                    "Using layer '{.val ",
                    layer_name,
                    "}' for count matrix",
                    message_type="success",
                )
                break
            except Exception as e:
                log_message(
                    "Failed to access layer '{.val ",
                    layer_name,
                    "}': {.val ",
                    str(e),
                    "}",
                    message_type="warning",
                )

    # If no specific layer found, use first available layer
    if count_matrix is None:
        layer_name = list(adata.layers.keys())[0]
        log_message(
            "Using first available layer '{.val ",
            layer_name,
            "}' for count matrix...",
            message_type="info",
        )
        count_matrix = adata.layers[layer_name]
        gene_names = adata.var_names
        var_df = adata.var.copy()

# Final check
if count_matrix is None:
    log_message(
        "Error: Could not extract count matrix from any source", message_type="error"
    )
    sys.exit(1)

# Check and convert matrix format for R compatibility
log_message("Checking matrix format...", message_type="info")
if sparse.issparse(count_matrix):
    log_message("Count matrix is sparse, format: {.val ", count_matrix.format, "}")
    # Convert to CSC format for better R compatibility
    if count_matrix.format != "csc":
        log_message(
            "Converting to CSC format for R compatibility...", message_type="running"
        )
        count_matrix = count_matrix.tocsc()
else:
    log_message(
        "Count matrix is dense, converting to sparse CSC format...",
        message_type="running",
    )
    count_matrix = sparse.csc_matrix(count_matrix)

log_message("Count matrix format: {.val ", count_matrix.format, "}")

# Create a new AnnData object with count matrix in X
# This ensures X is in the correct format for scop's adata_to_srt
log_message("Creating compatible AnnData object...", message_type="running")

adata_new = sc.AnnData(X=count_matrix, obs=adata.obs.copy(), var=var_df.copy())
adata_new.var_names = gene_names
adata_new.obs_names = adata.obs_names

# Save compatible h5ad file
compatible_h5ad_path = os.path.join(res_dir, "BTS_atlas_compatible.h5ad")
log_message("Saving compatible h5ad file to {.path ", compatible_h5ad_path, "}")
adata_new.write(compatible_h5ad_path, compression="gzip")

log_message("Compatible h5ad file saved successfully", message_type="success")
log_message(
    "You can now use this file with scop's adata_to_srt function", message_type="info"
)

# Also save metadata as CSV for R script
log_message("Saving metadata as CSV...", message_type="running")

# Extract metadata
metadata = adata.obs.copy()
metadata.index.name = "cell_barcode"

log_message("Metadata shape: {.val ", metadata.shape, "}")
log_message("Metadata columns: {.val ", list(metadata.columns), "}")

# Save metadata as CSV
metadata_csv_path = os.path.join(res_dir, "metadata.csv")
log_message("Saving metadata to {.path ", metadata_csv_path, "}")
metadata.to_csv(metadata_csv_path)

log_message("Metadata saved successfully", message_type="success")

# Save gene information if available
if var_df.shape[1] > 0:
    gene_info = var_df.copy()
    gene_info.index.name = "gene"
    gene_info_csv_path = os.path.join(res_dir, "gene_info.csv")
    log_message("Saving gene information to {.path ", gene_info_csv_path, "}")
    gene_info.to_csv(gene_info_csv_path)
    log_message("Gene information saved successfully", message_type="success")

log_message("Data conversion completed!", message_type="success")
log_message(
    "Compatible h5ad file: {.path ", compatible_h5ad_path, "}", message_type="info"
)
