source("code/functions/hic.R")

results_dir <- "results/hi-c/"
figures <- "results/figures/"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures, recursive = TRUE, showWarnings = FALSE)

color_palette_method <- c(
  "GENIE3" = "#cc3366",
  "LEAP" = "#ff9900",
  "PPCOR" = "#339999",
  "HARNexus" = "#005ea3"
)

log_message("Loading data...")
if (!file.exists(paste0(results_dir, "/data/input_data.Rdata"))) {
  load("data/astro_CBC_1kb_geneupdown_2kb.Rdata")
  load("data/Astro_network_table_list.Rdata")

  cbc_astro_hic <- overlap_info[, c(2, 6, 3, 9)]
  colnames(cbc_astro_hic) <- c("HAR", "gene_hic", "score_hic", "overlap")
  cbc_astro_hic <- na.omit(cbc_astro_hic)

  tfdb_result <- read.delim(
    "data/AnimalTFDB_result.txt",
    stringsAsFactors = FALSE
  )
  tfdb_result <- tfdb_result[, c(1, 3, 7)]
  colnames(tfdb_result) <- c("TF", "HAR", "score_TFBS")

  network_table <- list(
    "GENIE3" = network_table_genie3,
    "LEAP" = network_table_leap,
    "PPCOR" = network_table_ppcor,
    "HARNexus" = network_table_infercsn
  )

  save(
    cbc_astro_hic,
    tfdb_result,
    network_table,
    network_table_genie3,
    network_table_leap,
    network_table_ppcor,
    network_table_infercsn,
    file = paste0(results_dir, "/data/input_data.Rdata")
  )
} else {
  load(paste0(results_dir, "/data/input_data.Rdata"))
}

if (!file.exists(paste0(results_dir, "/data/intersection_results.Rdata"))) {
  intersection_results <- calculate_intersection_trend(
    network_table,
    tfdb_result,
    cbc_astro_hic,
    max_edges = 9000,
    step_size = 500
  )
  save(
    intersection_results,
    file = paste0(results_dir, "/data/intersection_results.Rdata")
  )
} else {
  load(paste0(results_dir, "/data/intersection_results.Rdata"))
}

coverage_plot <- coverage_trend_plot(
  network_table,
  color_palette_method
)
ggsave(
  paste0(figures, "/fig.s2b_intersection_ratio_trend.pdf"),
  coverage_plot + theme(
    legend.position = "right",
  ),
  width = 6,
  height = 3
)

intersection_ratio_plot <- intersection_ratio_trend_plot(
  intersection_results,
  color_palette_method
)
ggsave(
  paste0(figures, "/fig.2b_intersection_ratio_trend.pdf"),
  intersection_ratio_plot,
  width = 5,
  height = 2.5
)

combined_plot <- coverage_plot +
  intersection_ratio_plot +
  plot_annotation(
    tag_levels = "a"
  ) +
  plot_layout(
    ncol = 2,
    guides = "collect"
  ) & theme(
  legend.position = "bottom"
)

ggsave(
  paste0(figures, "/fig.S2_intersection_ratio.pdf"),
  combined_plot,
  width = 7.5,
  height = 3.5
)

combined_plot_2 <- coverage_plot +
  intersection_ratio_plot +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  ) +
  plot_layout(
    ncol = 2,
    guides = "collect"
  ) & theme(
  legend.position = "bottom"
)

ggsave(
  paste0(figures, "/fig.S2_intersection_ratio_2.pdf"),
  combined_plot_2,
  width = 7.5, height = 3.5
)


weight_plots_list <- create_weight_plots(
  network_table,
  intersection_results,
  thresholds = 1000
)
weight_plots_list[["GENIE3_1000"]] <- weight_plots_list[["GENIE3_1000"]] +
  labs(
    x = "Edges"
  )
weight_plots_list[["LEAP_1000"]] <- weight_plots_list[["LEAP_1000"]] +
  labs(
    y = "",
    x = "Edges"
  )
weight_plots_list[["PPCOR_1000"]] <- weight_plots_list[["PPCOR_1000"]] +
  labs(
    y = "",
    x = "Edges"
  )
weight_plots_list[["HARNexus"]] <- weight_plots_list[["HARNexus"]] +
  labs(
    y = "",
    x = "Edges"
  )
weight_plots <- weight_plots_list[["GENIE3_1000"]] +
  weight_plots_list[["LEAP_1000"]] +
  weight_plots_list[["PPCOR_1000"]] +
  weight_plots_list[["HARNexus"]] +
  plot_annotation(
    tag_levels = "a"
  ) +
  plot_layout(
    ncol = 4,
    guides = "collect"
  ) &
  theme(legend.position = "bottom")


comparison_plot <- create_comparison_barplot(
  intersection_results,
  thresholds = 1000,
  color_palette_method
)

comparison_plot <- comparison_plot +
  theme(
    legend.position = "bottom"
  )

venn_plots_list <- create_venn_diagrams(
  intersection_results,
  thresholds = 1000,
  color_palette_method
)

venn_plot <- venn_plots_list[[1]]
venn_plot <- venn_plot +
  theme(
    legend.position = "none"
  )
top_row <- comparison_plot + venn_plot +
  plot_layout(
    ncol = 2,
    widths = c(0.45, 0.55)
  )

final_plot <- top_row / weight_plots +
  plot_layout(
    heights = c(0.65, 0.35)
  ) +
  plot_annotation(
    tag_levels = "a"
  ) &
  theme(
    text = element_text(size = 15)
  )

ggsave(
  paste0(figures, "/fig.2b_network_comparison.pdf"),
  final_plot,
  width = 12,
  height = 10
)

ggsave(
  paste0(figures, "/fig.2b_network_comparison.png"),
  final_plot & theme(text = element_text(family = "Times New Roman")),
  width = 12,
  height = 10,
  dpi = 600,
  bg = "white"
)

final_plot_2 <- top_row / weight_plots +
  plot_layout(
    heights = c(0.65, 0.35)
  ) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  theme(
    text = element_text(size = 15)
  )

ggsave(
  paste0(figures, "/fig.2b_network_comparison_2.pdf"),
  final_plot_2,
  width = 12,
  height = 10
)
ggsave(
  paste0(figures, "/fig.2b_network_comparison_2.png"),
  final_plot_2 & theme(text = element_text(family = "Times New Roman")),
  width = 12,
  height = 10,
  dpi = 600,
  bg = "white"
)
