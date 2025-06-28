source("code/functions/packages.R")

result <- readRDS("data/networks/reorganized_networks.rds")

network_stats <- function(network_data) {
  stats <- data.frame()

  for (stage in names(network_data)) {
    for (region in names(network_data[[stage]])) {
      for (cell_type in names(network_data[[stage]][[region]])) {
        network <- network_data[[stage]][[region]][[cell_type]]
        if (!is.null(network) && nrow(network) > 0) {
          edge_count <- nrow(network)
          unique_tfs <- length(unique(network$regulator))
          unique_targets <- length(unique(network$target))

          stats <- rbind(
            stats, data.frame(
              Stage = stage,
              Region = region,
              Cell_type = cell_type,
              "Edges_count" = edge_count,
              "TFs_count" = unique_tfs,
              "Genes_count" = unique_targets
            )
          )
        }
      }
    }
  }

  stats$Stage <- factor(stats$Stage,
    levels = paste0("S", 1:10),
    ordered = TRUE
  )

  return(stats)
}


create_single_heatmap <- function(
    data,
    cell_type,
    metric = "Edges_count",
    legend_label = NULL) {
  legend_labels <- list(
    "Edges_count" = "Edge counts",
    "TFs_count" = "TF counts",
    "Genes_count" = "Gene counts"
  )

  p <- ggplot(
    data,
    aes(
      x = Region,
      y = Stage,
      fill = .data[[metric]]
    )
  ) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_viridis(
      option = "viridis",
      na.value = "white"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(angle = 0),
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "lightgrey", color = "white")
    ) +
    labs(
      title = cell_type,
      x = "Brain Region",
      y = "Developmental Stage",
      fill = legend_labels[[metric]]
    )

  return(p)
}

create_heatmap <- function(stats, metric = "Edges_count") {
  cell_types <- sort(unique(stats$Cell_type))

  value_range <- range(stats[[metric]], na.rm = TRUE)

  plot_list <- list()
  for (cell_type in cell_types) {
    cell_data <- stats[stats$Cell_type == cell_type, ]
    p <- create_single_heatmap(cell_data, cell_type, metric)
    # Add consistent limits to scale_fill_viridis
    p <- p + scale_fill_viridis(
      option = "viridis",
      na.value = "white",
      limits = value_range
    )
    plot_list[[cell_type]] <- p
  }

  combined_plot <- wrap_plots(plot_list, ncol = 3)

  return(combined_plot)
}

dir.create("results/networks", recursive = TRUE, showWarnings = FALSE)

network_stats <- network_stats(result)
head(network_stats)

write.xlsx(
  network_stats,
  file = "results/networks/network_statistics.xlsx",
  sheetName = "Network Stats",
  rowNames = FALSE,
  colNames = TRUE,
  borders = "surrounding"
)

metrics <- c("Edges_count", "TFs_count", "Genes_count")
for (metric in metrics) {
  log_message("Processing ", metric, "...")
  combined_plot <- create_heatmap(network_stats, metric)
  if (metric == "Edges_count") {
    file_name <- "fig.4a"
  } else if (metric == "TFs_count") {
    file_name <- "fig.S6"
  } else if (metric == "Genes_count") {
    file_name <- "fig.S7"
  }
  ggsave(
    paste0("results/networks/", file_name, ".pdf"),
    combined_plot,
    width = 13,
    height = 8
  )
}
