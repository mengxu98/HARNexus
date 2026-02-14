#!/bin/bash

set -e

source "code/functions/utils.sh"

log_message "Running all plotting scripts"
echo ""

plotting_dir="code/plotting"
run_r_script "$plotting_dir" "har_tf_pairs_plots.R" "HAR-TF pairs analysis visualization"
run_r_script "$plotting_dir" "tfs_enrichment_analysis.R" "TF enrichment analysis visualization"
run_r_script "$plotting_dir" "network_comparison.R" "Network method comparison"
run_r_script "$plotting_dir" "GSE97942_plotting.R" "GSE97942 dataset visualization"
run_r_script "$plotting_dir" "hic_intersect_analysis.R" "Hi-C intersection analysis"
run_r_script "$plotting_dir" "motifs_processing_plots.R" "Motif processing results visualization"
run_r_script "$plotting_dir" "aucell.R" "AUCell analysis visualization"
run_r_script "$plotting_dir" "datasets.R" "Dataset integration visualization"
run_r_script "$plotting_dir" "lisi_plot.R" "LISI visualization"
run_r_script "$plotting_dir" "altas_statistics.R" "Atlas statistics visualization"
run_python_script "$plotting_dir" "pfc_celltype_upset.py" "PFC cell type UpSet plot"
run_r_script "$plotting_dir" "pfc_astrocytes_analysis.R" "PFC astrocytes analysis visualization"
run_r_script "$plotting_dir" "similarity_plot.R" "Similarity analysis visualization"
run_python_script "$plotting_dir" "network_sankey.py" "Network Sankey diagram (static)"
run_python_script "$plotting_dir" "network_sankey_interactive.py" "Network Sankey diagram (interactive)"
run_r_script "$plotting_dir" "gse192774_evolution_plots.R" "Comprehensive evolution analysis visualization"
run_r_script "$plotting_dir" "gse192774_network_visualize.R" "Network visualization (GSE192774)"
run_r_script "$plotting_dir" "gse192774_atac_daccre_plots.R" "DAcCRE visualization"
run_r_script "$plotting_dir" "gse192774_similarity_heatmap.R" "Network similarity heatmap"
run_r_script "$plotting_dir" "gse192774_dim_plot.R" "Dimensionality reduction plots"
run_r_script "$plotting_dir" "gse192774_scatter.R" "Scatter plots"

log_success "All plotting scripts completed!"
echo ""
