#!/bin/bash

set -e

source "code/functions/utils.sh"

log_message "Running all plotting scripts"
echo ""

plotting_dir="code/plotting"

# Figure 1a
run_r_script "$plotting_dir" "har_tf_pairs_plots.R" "HAR-TF pairs analysis visualization"
# Figure 1b
run_r_script "$plotting_dir" "tfs_enrichment_analysis.R" "TF enrichment analysis visualization"
# Figure 1c-e
run_r_script "$plotting_dir" "motifs_processing_plots.R" "Motif processing results visualization"

# Supplementary figure 1a-e
run_r_script "$plotting_dir" "GSE97942_plotting.R" "GSE97942 dataset visualization"
# Figure 2a-c, Supplementary figure 1f-g
run_r_script "$plotting_dir" "network_comparison.R" "Network method comparison"
# Figure 2b, 2d-e
run_r_script "$plotting_dir" "hic_intersect_analysis.R" "Hi-C intersection analysis"

# Figure 3, Supplementary figure 2a-b
run_r_script "$plotting_dir" "aucell.R" "AUCell analysis visualization"

# Figure 4a-c, Supplementary figure 3a, 3c-f, Supplementary figure 4a-b
run_r_script "$plotting_dir" "datasets.R" "Dataset integration visualization"
# Supplementary figure 3b
run_r_script "$plotting_dir" "lisi_plot.R" "LISI visualization"

# Figure 4d, Supplementary figure 5a-b
run_r_script "$plotting_dir" "altas_statistics.R" "Atlas statistics visualization"
# Figure 4e, Figure 6g
run_python_script "$plotting_dir" "network_sankey.py" "Network Sankey diagram (static)"
# run_python_script "$plotting_dir" "network_sankey_interactive.py" "Network Sankey diagram (interactive)"
# Figure 4f-g, Supplementary figure 6a-b
run_r_script "$plotting_dir" "similarity_plot.R" "Similarity analysis visualization"
# Supplementary figure 6c-d
run_python_script "$plotting_dir" "pfc_celltype_upset.py" "PFC cell type UpSet plot"

# Figure 5a-b, Supplementary figure 7a-b
run_r_script "$plotting_dir" "gse192774_dim_plot.R" "Dimensionality reduction plots"
# Figure 5c, Supplementary figure 7c
run_r_script "$plotting_dir" "gse192774_scatter.R" "Scatter plots"
# Figure 5f-g, Supplementary figure 7f-g, Supplementary figure 8a-b
run_r_script "$plotting_dir" "gse192774_similarity_heatmap.R" "Network similarity heatmap"
# Figure 5h-i, Supplementary figure 7h-i
run_r_script "$plotting_dir" "gse192774_evolution_plots.R" "Comprehensive evolution analysis visualization"
# Figure 5j, Supplementary figure 7j
run_r_script "$plotting_dir" "gse192774_network_visualize.R" "Network visualization (GSE192774)"
# Not used in both of main and supplementary figures
run_r_script "$plotting_dir" "gse192774_atac_daccre_plots.R" "DAcCRE visualization"

# Figure 6a-f, Supplementary figure 9a-d
run_r_script "$plotting_dir" "pfc_astrocytes_analysis.R" "PFC astrocytes analysis visualization"

log_success "All plotting scripts running completed!"
echo ""
