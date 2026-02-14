source("code/functions/hic.R")
source("code/functions/prepare_env.R")

fig_dir <- check_dir("figures/hic")
edge_num <- 500

log_message("Loading data...")
network_list <- readRDS(
  "results/gse97942/astro_4methods_networks.rds"
)

intersection_results <- readRDS(
  "results/hic/intersection_results.rds"
)

hic_data <- read.csv(
  "results/hic/HAR_gene_HiC_supported.csv",
  header = TRUE
)

tfdb_data <- read.csv(
  "results/har_tf/human/har_tf_pairs_scores.csv"
)
tfdb_data <- tfdb_data[, c("TF", "Sequence", "har", "jaspar_score")]
colnames(tfdb_data) <- c("TF", "HAR", "HAR_ID", "Score")

intersection_ratio_plot <- intersection_ratio_trend_plot(
  intersection_results,
  color_methods
)
ggsave(
  paste0(fig_dir, "/intersection_ratio.pdf"),
  intersection_ratio_plot,
  width = 5.2,
  height = 2.5
)

weight_plots <- weight_ratio_plots(
  network_list = network_list,
  intersection_results = intersection_results,
  thresholds = edge_num,
  color_methods = color_methods,
  tfdb_data = tfdb_data,
  hic_data = hic_data,
  ratio_plot_type = "pie",
  top_n_validated = 0,
  ncol = 4,
  brain_pattern = brain_pattern,
  go_ontology = "ALL"
)
ggsave(
  paste0(fig_dir, "/weight_plots.pdf"),
  weight_plots,
  width = 10.5,
  height = 11
)

comparison_plot <- comparison_barplot(
  intersection_results,
  thresholds = edge_num,
  color_methods
)
ggsave(
  paste0(fig_dir, "/comparison_barplot.pdf"),
  comparison_plot,
  width = 5,
  height = 2.5
)

venn_plots_list <- create_venn_diagrams(
  intersection_results,
  type = "all",
  thresholds = edge_num,
  color_methods
)
venn_plot <- venn_plots_list[[1]]
ggsave(
  paste0(fig_dir, "/venn_plot.pdf"),
  venn_plot,
  width = 6,
  height = 4.5
)
