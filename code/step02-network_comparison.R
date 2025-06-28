source("code/functions/packages.R")
source("code/functions/network_comparison.R")

figures_paper <- "results/figures_paper"
figures_zhang <- "results/figures"

color_palette_method <- c(
  "GENIE3" = "#cc3366",
  "LEAP" = "#ff9900",
  "PPCOR" = "#339999",
  "HARNexus" = "#005ea3"
)

color_palette_edge <- c(
  "500 edges" = "#84bbee",
  "1000 edges" = "#4993d8",
  "2000 edges" = "#2a7bc7",
  "HARNexus" = "#005ea3"
)

methods_order <- c("GENIE3", "PPCOR", "LEAP", "HARNexus")

load("data/Astro_network_table_list.Rdata")

tfs_harnexus <- as.character(unique(network_table_infercsn$regulator))
genes_infercsn <- as.character(unique(network_table_infercsn$target))
tfs_genie3 <- as.character(unique(network_table_genie3$regulator))
genes_genie3 <- as.character(unique(network_table_genie3$target))
tfs_leap <- as.character(unique(network_table_leap$regulator))
genes_leap <- as.character(unique(network_table_leap$target))
tfs_ppcor <- as.character(unique(network_table_ppcor$regulator))
genes_ppcor <- as.character(unique(network_table_ppcor$target))

nodes <- reorder_data(
  c(
    length(c(tfs_genie3, genes_genie3)),
    length(c(tfs_ppcor, genes_ppcor)),
    length(c(tfs_leap, genes_leap)),
    length(c(tfs_harnexus, genes_infercsn))
  )
)
edges <- reorder_data(
  c(
    nrow(network_table_genie3),
    nrow(network_table_ppcor),
    nrow(network_table_leap),
    nrow(network_table_infercsn)
  )
)
tfs <- reorder_data(
  c(
    length(tfs_genie3),
    length(tfs_ppcor),
    length(tfs_leap),
    length(tfs_harnexus)
  )
)
genes <- reorder_data(
  c(
    length(genes_genie3),
    length(genes_ppcor),
    length(genes_leap),
    length(genes_infercsn)
  )
)

nodes_breaks <- c(200, 1500)
edges_breaks <- c(600, 60000)
tfs_breaks <- c(10, 50)
genes_breaks <- c(200, 1500)

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
  color_palette_method,
  nodes_breaks
)

p2 <- create_split_plot(
  edges_data,
  "Number of edges",
  "",
  color_palette_method,
  edges_breaks
)

p3 <- create_split_plot(
  tfs_data,
  "Number of TFs",
  "",
  color_palette_method,
  tfs_breaks
)

p4 <- create_split_plot(
  genes_data,
  "Number of genes",
  "",
  color_palette_method,
  genes_breaks
)

combined_plot <- p2 + p1 + p3 + p4 +
  plot_annotation(
    tag_levels = "a"
  ) +
  plot_layout(
    ncol = 4,
    guides = "collect"
  ) & theme(
  legend.position = "none"
)

ggsave(
  paste0(figures_paper, "/fig.S1_network_comparison.pdf"),
  combined_plot,
  width = 12,
  height = 3
)

ggsave(
  paste0(figures_paper, "/fig.S1_network_comparison.png"),
  combined_plot,
  width = 12,
  height = 3,
  dpi = 600
)
combined_plot_22 <- p2 + p1 + p3 + p4 +
  plot_annotation(
    title = "",
    theme = theme(
      plot.title = element_text(hjust = 0.5)
    ),
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  ) +
  plot_layout(
    ncol = 2,
    guides = "collect"
  ) & theme(
  legend.position = "none"
)
ggsave(
  paste0(figures_zhang, "/network_comparison-2x2.pdf"),
  combined_plot_22,
  width = 8,
  height = 6
)


k <- 1000
network_genie3_k <- network_table_genie3[1:1000, ]
network_ppcor_k <- network_table_ppcor[1:1000, ]
network_leap_k <- network_table_leap[1:1000, ]
network_infercsn_k <- network_table_infercsn

tfs_genie3_k <- as.character(unique(network_genie3_k$regulator))
tfs_ppcor_k <- as.character(unique(network_ppcor_k$regulator))
tfs_leap_k <- as.character(unique(network_leap_k$regulator))
tfs_infercsn_k <- as.character(unique(network_infercsn_k$regulator))

genes_genie3_k <- as.character(unique(network_genie3_k$target))
genes_ppcor_k <- as.character(unique(network_ppcor_k$target))
genes_leap_k <- as.character(unique(network_leap_k$target))
genes_infercsn_k <- as.character(unique(network_infercsn_k$target))

nodes_k <- reorder_data(
  c(
    length(c(tfs_genie3_k, genes_genie3_k)),
    length(c(tfs_ppcor_k, genes_ppcor_k)),
    length(c(tfs_leap_k, genes_leap_k)),
    length(c(tfs_infercsn_k, genes_infercsn_k))
  )
)

edges_k <- reorder_data(
  c(
    nrow(network_genie3_k),
    nrow(network_ppcor_k),
    nrow(network_leap_k),
    nrow(network_infercsn_k)
  )
)

tfs_k <- reorder_data(
  c(
    length(tfs_genie3_k),
    length(tfs_ppcor_k),
    length(tfs_leap_k),
    length(tfs_infercsn_k)
  )
)

genes_k <- reorder_data(
  c(
    length(genes_genie3_k),
    length(genes_ppcor_k),
    length(genes_leap_k),
    length(genes_infercsn_k)
  )
)

nodes_breaks_k <- c(150, 250)
edges_breaks_k <- c(400, 450)
tfs_breaks_k <- c(30, 45)
genes_breaks_k <- c(150, 250)

nodes_data_k <- create_split_data(
  methods_order,
  nodes_k,
  nodes_breaks_k
)

edges_data_k <- create_split_data(
  methods_order,
  edges_k,
  edges_breaks_k
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
  color_palette_method,
  nodes_breaks_k
)

edges_plot_k <- create_split_plot(
  edges_data_k,
  "Number of edges",
  "",
  color_palette_method,
  edges_breaks_k
)

tfs_plot_k <- create_split_plot(
  tfs_data_k,
  "Number of TFs",
  "",
  color_palette_method,
  tfs_breaks_k
)

genes_plot_k <- create_split_plot(
  genes_data_k,
  "Number of genes",
  "",
  color_palette_method,
  genes_breaks_k
)

plot_network_comparison <- edges_plot_k + nodes_plot_k +
  tfs_plot_k + genes_plot_k +
  plot_annotation(
    tag_levels = "a"
  ) +
  plot_layout(
    ncol = 4,
    guides = "collect"
  ) & theme(legend.position = "none")

ggsave(
  paste0(figures_paper, "/fig.2a_network_comparison_1000.pdf"),
  plot_network_comparison,
  width = 12,
  height = 3
)
