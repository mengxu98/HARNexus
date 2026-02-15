source("code/functions/hic.R")
source("code/functions/prepare_env.R")

res_dir <- check_dir("results/hic")

log_message("Loading data...")
network_list <- readRDS(
  "results/gse97942/astro_4methods_networks.rds"
)
plot_network_distribution(network_list$GENIE3)
plot_network_distribution(network_list$HARNexus)
plot_network_distribution(network_list$LEAP)
plot_network_distribution(network_list$PPCOR)

hic_data <- read.csv(
  "results/hic/HAR_gene_HiC_supported.csv",
  header = TRUE
)

tfdb_data <- read.csv(
  "results/har_tf/human/har_tf_pairs_scores.csv"
)
tfdb_data <- tfdb_data[, c("TF", "Sequence", "har", "jaspar_score")]
colnames(tfdb_data) <- c("TF", "HAR", "HAR_ID", "Score")

if (!file.exists(paste0(res_dir, "/intersection_results.rds"))) {
  intersection_results <- calculate_intersection_trend(
    network_list,
    tfdb_data = tfdb_data,
    hic_data = hic_data,
    max_edges = 10000,
    step_size = 500
  )
  saveRDS(
    intersection_results,
    paste0(res_dir, "/intersection_results.rds")
  )
}
