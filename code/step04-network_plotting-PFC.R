source("code/functions/utils.R")
source("code/functions/network_analysis.R")

results_dir <- check_dir("results/networks_plotting_pfc/")

if (!file.exists(file.path(results_dir, "pfc_network.rds"))) {
  log_message("Processing data...")
  network_data <- read.csv("data/networks/csv/network_data.csv")
  pfc_network <- network_data[network_data$Region == "PFC", ]
  saveRDS(
    pfc_network,
    file.path(
      results_dir, "pfc_network.rds"
    )
  )
} else {
  log_message("Loading data...")
  pfc_network <- readRDS(
    file.path(results_dir, "pfc_network.rds")
  )
}

stages <- sort(unique(pfc_network$Stage))
n_stages <- length(stages)
log_message(
  sprintf(
    "Found %d stages: %s",
    n_stages,
    paste(stages, collapse = ", ")
  )
)

log_message("Finding genes common to S8 and S9 but not in other stages...")

for (stage in stages) {
  stage_data <- pfc_network[pfc_network$Stage == stage, ]
  log_message(
    sprintf(
      "Stage %s: %d unique targets", stage, length(unique(stage_data$Target))
    )
  )
}

s8_s9_common_genes <- unique(
  unique(pfc_network[pfc_network$Stage == "S8", ]$Target),
  unique(pfc_network[pfc_network$Stage == "S9", ]$Target)
)
other_stages_genes <- unique(
  Reduce(
    union, lapply(
      setdiff(stages, c("S8", "S9")),
      function(stage) {
        unique(pfc_network[pfc_network$Stage == stage, ]$Target)
      }
    )
  )
)
s8_s9_unique_genes <- unique(setdiff(s8_s9_common_genes, other_stages_genes))

log_message(
  sprintf(
    "Found %d genes common to S8 and S9, %d unique to S8 and S9",
    length(s8_s9_common_genes),
    length(s8_s9_unique_genes)
  )
)

tf_lists <- lapply(
  stages, function(stage) {
    stage_data <- pfc_network[pfc_network$Stage == stage, ]
    unique(stage_data$TF)
  }
)
names(tf_lists) <- stages

target_lists <- lapply(
  stages, function(stage) {
    stage_data <- pfc_network[pfc_network$Stage == stage, ]
    unique(stage_data$Target)
  }
)
names(target_lists) <- stages

intersect_tfs <- Reduce(intersect, tf_lists)
intersect_targets <- Reduce(intersect, target_lists)

setdiff_targets_list <- lapply(
  stages, function(stage) {
    stage_data <- pfc_network[pfc_network$Stage == stage, ]
    stage_data_other <- pfc_network[pfc_network$Stage != stage, ]
    setdiff(stage_data$Target, stage_data_other$Target)
  }
)
names(setdiff_targets_list) <- stages
setdiff_targets <- Reduce(union, setdiff_targets_list)

all_targets <- Reduce(union, target_lists)

log_message(
  sprintf(
    "Found %d intersect TFs, and %d intersect targets...",
    length(intersect_tfs),
    length(intersect_targets)
  )
)

pfc_network_intersect <- pfc_network[pfc_network$Target %in% intersect_targets, ]

stages <- sort(unique(pfc_network$Stage))
n_stages <- length(stages)
log_message(
  sprintf(
    "Found %d stages: %s",
    n_stages,
    paste(stages, collapse = ", ")
  )
)

log_message("Plotting networks...")

plot_network_stages <- list()
for (stage in stages) {
  stage_data <- pfc_network_intersect[pfc_network_intersect$Stage == stage, ]

  p_combined <- plot_networks(
    stage_data,
    stage,
    max_edges_per_type = 200,
    split = FALSE,
    tfs = intersect_tfs,
    targets = intersect_targets,
    label_size = 5,
    show_tf_labels = TRUE,
    show_gene_labels = FALSE,
    node_colors = c(TF = "black", Gene = "gray50"),
    node_size = c(TF = 5, Gene = 3.5),
    node_color_celltypes = TRUE
  )
  plot_network_stages[[stage]] <- p_combined
}

combined_plot <- wrap_plots(
  plot_network_stages,
  ncol = 3,
  guides = "collect"
) & theme(legend.position = "right")

ggsave(
  file.path(results_dir, "fig.S8-networks_all_stages.pdf"),
  combined_plot,
  width = 20,
  height = 20
)

for (stage in stages) {
  stage_data <- pfc_network_intersect[pfc_network_intersect$Stage == stage, ]

  p_split <- plot_networks(
    stage_data,
    stage,
    max_edges_per_type = 500,
    split = TRUE,
    tfs = intersect_tfs,
    targets = intersect_targets
  )
  if (!is.null(p_split)) {
    if (stage == "S3") {
      width <- 10
      height <- 6
    } else {
      width <- 12
      height <- 15
    }
    ggsave(
      file.path(results_dir, sprintf("networks_split_%s.pdf", stage)),
      p_split,
      width = width,
      height = height
    )
  }
}
