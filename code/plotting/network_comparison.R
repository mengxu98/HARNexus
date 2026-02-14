source("code/functions/network_comparison.R")
source("code/functions/prepare_env.R")

k <- 500
exclude_tfs_in_genes <- TRUE

fig_dir <- check_dir("figures/gse97942/")

methods_order <- c("GENIE3", "HARNexus", "LEAP", "PPCOR")

network_list <- readRDS(
  "results/networks/gse97942/astro_4methods_networks.rds"
)

coverage_plot <- coverage_trend_plot(
  network_list,
  color_methods,
  text_size = 3.5
)
ggsave(
  paste0(fig_dir, "/intersection_ratio_trend.pdf"),
  coverage_plot,
  width = 4.5,
  height = 3.5
)

network_genie3 <- network_list$GENIE3
if (exclude_tfs_in_genes) {
  network_genie3 <- network_genie3[!network_genie3$target %in% network_genie3$regulator, ]
}
network_harnexus <- network_list$HARNexus
if (exclude_tfs_in_genes) {
  network_harnexus <- network_harnexus[!network_harnexus$target %in% network_harnexus$regulator, ]
}
network_ppcor <- network_list$PPCOR
if (exclude_tfs_in_genes) {
  network_ppcor <- network_ppcor[!network_ppcor$target %in% network_ppcor$regulator, ]
}
network_leap <- network_list$LEAP
if (exclude_tfs_in_genes) {
  network_leap <- network_leap[!network_leap$target %in% network_leap$regulator, ]
}

tfs_harnexus <- as.character(unique(network_harnexus$regulator))
tfs_genie3 <- as.character(unique(network_genie3$regulator))
tfs_leap <- as.character(unique(network_leap$regulator))
tfs_ppcor <- as.character(unique(network_ppcor$regulator))

genes_genie3 <- as.character(unique(network_genie3$target))
genes_harnexus <- as.character(unique(network_harnexus$target))
genes_leap <- as.character(unique(network_leap$target))
genes_ppcor <- as.character(unique(network_ppcor$target))

nodes <- reorder_data(
  c(
    length(c(tfs_genie3, genes_genie3)),
    length(c(tfs_harnexus, genes_harnexus)),
    length(c(tfs_ppcor, genes_ppcor)),
    length(c(tfs_leap, genes_leap))
  )
)
edges <- reorder_data(
  c(
    nrow(network_genie3),
    nrow(network_harnexus),
    nrow(network_leap),
    nrow(network_ppcor)
  )
)
tfs <- reorder_data(
  c(
    length(tfs_genie3),
    length(tfs_harnexus),
    length(tfs_leap),
    length(tfs_ppcor)
  )
)
genes <- reorder_data(
  c(
    length(genes_genie3),
    length(genes_harnexus),
    length(genes_leap),
    length(genes_ppcor)
  )
)

nodes_breaks <- NULL
edges_breaks <- NULL
tfs_breaks <- NULL
genes_breaks <- NULL

nodes_data <- create_split_data(
  methods_order,
  nodes,
  nodes_breaks
)
edges_data <- create_split_data(
  methods_order,
  edges,
  edges_breaks
)
tfs_data <- create_split_data(
  methods_order,
  tfs,
  tfs_breaks
)
genes_data <- create_split_data(
  methods_order,
  genes,
  genes_breaks
)

p1 <- create_split_plot(
  nodes_data,
  "Number of nodes",
  "",
  color_methods,
  nodes_breaks
)

p2 <- create_split_plot(
  edges_data,
  "Number of edges",
  "",
  color_methods,
  edges_breaks
)

p3 <- create_split_plot(
  tfs_data,
  "Number of TFs",
  "",
  color_methods,
  tfs_breaks
)

p4 <- create_split_plot(
  genes_data,
  "Number of genes",
  "",
  color_methods,
  genes_breaks
)

combined_plot <- p2 + p1 + p3 + p4 +
  plot_layout(
    ncol = 2,
    guides = "collect"
  )

ggsave(
  file.path(fig_dir, "network_comparison.pdf"),
  combined_plot,
  width = 6.2,
  height = 5
)


network_genie3_k <- network_genie3[1:k, ]
network_harnexus_k <- network_harnexus[1:k, ]
network_ppcor_k <- network_ppcor[1:k, ]
network_leap_k <- network_leap[1:k, ]

tfs_genie3_k <- as.character(unique(network_genie3_k$regulator))
tfs_harnexus_k <- as.character(unique(network_harnexus_k$regulator))
tfs_ppcor_k <- as.character(unique(network_ppcor_k$regulator))
tfs_leap_k <- as.character(unique(network_leap_k$regulator))

genes_genie3_k <- as.character(unique(network_genie3_k$target))
genes_harnexus_k <- as.character(unique(network_harnexus_k$target))
genes_ppcor_k <- as.character(unique(network_ppcor_k$target))
genes_leap_k <- as.character(unique(network_leap_k$target))

nodes_k <- reorder_data(
  c(
    length(c(tfs_genie3_k, genes_genie3_k)),
    length(c(tfs_harnexus_k, genes_harnexus_k)),
    length(c(tfs_leap_k, genes_leap_k)),
    length(c(tfs_ppcor_k, genes_ppcor_k))
  )
)

tfs_k <- reorder_data(
  c(
    length(tfs_genie3_k),
    length(tfs_harnexus_k),
    length(tfs_leap_k),
    length(tfs_ppcor_k)
  )
)

genes_k <- reorder_data(
  c(
    length(genes_genie3_k),
    length(genes_harnexus_k),
    length(genes_leap_k),
    length(genes_ppcor_k)
  )
)

nodes_k
tfs_k
genes_k
nodes_breaks_k <- NULL
tfs_breaks_k <- NULL
genes_breaks_k <- NULL

nodes_data_k <- create_split_data(
  methods_order,
  nodes_k,
  nodes_breaks_k
)

tfs_data_k <- create_split_data(
  methods_order,
  tfs_k,
  tfs_breaks_k
)

genes_data_k <- create_split_data(
  methods_order,
  genes_k,
  genes_breaks_k
)

nodes_plot_k <- create_split_plot(
  nodes_data_k,
  "Number of nodes",
  "",
  color_methods,
  nodes_breaks_k
)

tfs_plot_k <- create_split_plot(
  tfs_data_k,
  "Number of TFs",
  "",
  color_methods,
  tfs_breaks_k
)

genes_plot_k <- create_split_plot(
  genes_data_k,
  "Number of genes",
  "",
  color_methods,
  genes_breaks_k
)

plot_network_comparison <- nodes_plot_k +
  tfs_plot_k + genes_plot_k +
  plot_annotation(
    tag_levels = "A"
  ) +
  plot_layout(
    guides = "collect"
  )

ggsave(
  file.path(fig_dir, paste0("network_comparison_", k, ".pdf")),
  plot_network_comparison,
  width = 8,
  height = 2.5
)
