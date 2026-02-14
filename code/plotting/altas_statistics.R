source("code/functions/prepare_env.R")

log_message("Loading network statistics data...")
network_stats <- read.csv(
  "results/networks/analysis/network_statistics.csv",
  stringsAsFactors = FALSE
)

create_single_heatmap <- function(
    data,
    cell_type,
    all_regions,
    all_stages,
    metric = "Edges_count",
    show_x_axis = TRUE,
    show_y_axis = TRUE,
    legend_label = NULL) {
  legend_labels <- list(
    "Edges_count" = "Edge counts",
    "TFs_count" = "TF counts",
    "Genes_count" = "Gene counts"
  )
  complete_grid <- expand.grid(
    Region = all_regions,
    Stage = all_stages,
    stringsAsFactors = FALSE
  )

  if (nrow(data) > 0) {
    data_subset <- data[, c("Region", "Stage", metric)]
    plot_data <- complete_grid %>%
      left_join(data_subset, by = c("Region", "Stage"))
    plot_data[[metric]][is.na(plot_data[[metric]])] <- 0
  } else {
    plot_data <- complete_grid
    plot_data[[metric]] <- 0
  }

  plot_data_for_plot <- plot_data
  plot_data_for_plot[[metric]][plot_data_for_plot[[metric]] == 0] <- NA

  plot_data_for_plot$Stage <- factor(
    plot_data_for_plot$Stage,
    levels = all_stages,
    ordered = TRUE
  )

  p <- ggplot(
    plot_data_for_plot,
    aes(
      x = Region,
      y = Stage,
      fill = .data[[metric]]
    )
  ) +
    geom_tile(color = "grey80", linewidth = 0.2) +
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(
      axis.text.x = if (show_x_axis) {
        element_text(angle = 60, hjust = 1)
      } else {
        element_blank()
      },
      axis.text.y = if (show_y_axis) {
        element_text(angle = 0)
      } else {
        element_blank()
      },
      axis.title.x = if (show_x_axis) {
        element_text()
      } else {
        element_blank()
      },
      axis.title.y = if (show_y_axis) {
        element_text()
      } else {
        element_blank()
      },
      axis.ticks.x = if (show_x_axis) {
        element_line()
      } else {
        element_blank()
      },
      axis.ticks.y = if (show_y_axis) {
        element_line()
      } else {
        element_blank()
      },
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(
        hjust = 0,
        color = if (cell_type %in% names(color_celltypes)) {
          color_celltypes[[cell_type]]
        } else {
          "black"
        }
      ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "lightgrey", color = "white"),
      plot.margin = margin(1, 1, 1, 1, "pt")
    ) +
    labs(
      subtitle = cell_type,
      x = if (show_x_axis) "Brain region" else NULL,
      y = if (show_y_axis) "Stage" else NULL,
      fill = legend_labels[[metric]]
    )

  return(p)
}

create_single_heatmap_region_stage <- function(
    data,
    cell_type,
    all_regions,
    all_stages,
    metric = "Edges_count",
    show_x_axis = TRUE,
    show_y_axis = TRUE) {
  legend_labels <- list(
    "Edges_count" = "Edge counts",
    "TFs_count" = "TF counts",
    "Genes_count" = "Gene counts"
  )
  complete_grid <- expand.grid(
    Region = all_regions,
    Stage = all_stages,
    stringsAsFactors = FALSE
  )

  if (nrow(data) > 0) {
    data_subset <- data[, c("Region", "Stage", metric)]
    plot_data <- complete_grid %>%
      left_join(data_subset, by = c("Region", "Stage"))
    plot_data[[metric]][is.na(plot_data[[metric]])] <- 0
  } else {
    plot_data <- complete_grid
    plot_data[[metric]] <- 0
  }

  plot_data_for_plot <- plot_data
  plot_data_for_plot[[metric]][plot_data_for_plot[[metric]] == 0] <- NA

  plot_data_for_plot$Stage <- factor(
    plot_data_for_plot$Stage,
    levels = all_stages,
    ordered = TRUE
  )
  plot_data_for_plot$Region <- factor(
    plot_data_for_plot$Region,
    levels = rev(all_regions),
    ordered = TRUE
  )

  p <- ggplot(
    plot_data_for_plot,
    aes(
      x = Stage,
      y = Region,
      fill = .data[[metric]]
    )
  ) +
    geom_tile(color = "grey80", linewidth = 0.2) +
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(
      axis.text.x = if (show_x_axis) {
        element_text(angle = 60, hjust = 1)
      } else {
        element_blank()
      },
      axis.text.y = if (show_y_axis) {
        element_text(angle = 0, size = rel(1))
      } else {
        element_blank()
      },
      axis.title.x = if (show_x_axis) element_text() else element_blank(),
      axis.title.y = if (show_y_axis) element_text() else element_blank(),
      axis.ticks.x = if (show_x_axis) element_line() else element_blank(),
      axis.ticks.y = if (show_y_axis) element_line() else element_blank(),
      plot.subtitle = element_text(
        hjust = 0.5,
        color = if (cell_type %in% names(color_celltypes)) {
          color_celltypes[[cell_type]]
        } else {
          "black"
        }
      ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(1, 1, 1, 1, "pt")
    ) +
    labs(
      subtitle = cell_type,
      x = if (show_x_axis) "Stage" else NULL,
      y = if (show_y_axis) "Brain region" else NULL,
      fill = legend_labels[[metric]]
    )

  return(p)
}

create_heatmap <- function(stats, metric = "Edges_count") {
  cell_types <- sort(unique(stats$Cell_type))

  all_regions <- sort(unique(stats$Region))
  all_stages <- sort(unique(stats$Stage))

  stage_numbers <- as.numeric(gsub("S", "", all_stages))
  all_stages <- all_stages[order(stage_numbers)]

  non_zero_values <- stats[[metric]][stats[[metric]] > 0]
  if (length(non_zero_values) > 0) {
    color_range <- range(non_zero_values, na.rm = TRUE)
  } else {
    color_range <- range(stats[[metric]], na.rm = TRUE)
  }

  plot_list <- list()
  n_cols <- 3
  n_rows <- ceiling(length(cell_types) / n_cols)

  for (i in seq_along(cell_types)) {
    cell_type <- cell_types[i]
    cell_data <- stats[stats$Cell_type == cell_type, ]

    show_y_axis <- ((i - 1) %% n_cols) == 0
    show_x_axis <- i > (n_rows - 1) * n_cols

    p <- create_single_heatmap(
      cell_data,
      cell_type,
      all_regions,
      all_stages,
      metric,
      show_x_axis = show_x_axis,
      show_y_axis = show_y_axis
    )

    p <- p + scale_fill_viridis(
      option = "viridis",
      na.value = "white",
      limits = color_range,
      oob = scales::oob_squish,
      labels = scales::number
    )
    plot_list[[cell_type]] <- p
  }

  combined_plot <- wrap_plots(plot_list, ncol = 3) +
    plot_layout(guides = "collect") &
    theme(
      legend.position = "bottom",
      legend.key.size = unit(0.28, "cm"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7, angle = 60, hjust = 1),
      plot.margin = margin(1, 1, 1, 1, "pt")
    )

  return(combined_plot)
}

create_heatmap_wide_region_stage <- function(stats, metric = "Edges_count") {
  cell_types <- sort(unique(stats$Cell_type))
  all_regions <- sort(unique(stats$Region))
  all_stages <- sort(unique(stats$Stage))
  stage_numbers <- as.numeric(gsub("S", "", all_stages))
  all_stages <- all_stages[order(stage_numbers)]

  non_zero_values <- stats[[metric]][stats[[metric]] > 0]
  if (length(non_zero_values) > 0) {
    color_range <- range(non_zero_values, na.rm = TRUE)
  } else {
    color_range <- range(stats[[metric]], na.rm = TRUE)
  }

  n_celltypes <- length(cell_types)
  n_cols <- min(9L, n_celltypes)

  plot_list <- list()
  for (i in seq_along(cell_types)) {
    cell_type <- cell_types[i]
    cell_data <- stats[stats$Cell_type == cell_type, ]
    show_y_axis <- (i - 1) %% n_cols == 0
    show_x_axis <- i > n_celltypes - n_cols
    p <- create_single_heatmap_region_stage(
      cell_data,
      cell_type,
      all_regions,
      all_stages,
      metric,
      show_x_axis = show_x_axis,
      show_y_axis = show_y_axis
    )
    p <- p + scale_fill_viridis(
      option = "viridis",
      na.value = "white",
      limits = color_range,
      oob = scales::oob_squish,
      labels = scales::number
    )
    plot_list[[cell_type]] <- p
  }

  combined_plot <- wrap_plots(plot_list, ncol = n_cols) +
    plot_layout(guides = "collect") &
    theme(
      legend.position = "bottom",
      legend.key.size = unit(0.28, "cm"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7, angle = 60, hjust = 1),
      plot.margin = margin(1, 1, 1, 1, "pt")
    )

  return(combined_plot)
}

fig_dir <- check_dir("figures/networks")

metrics <- c("Edges_count", "TFs_count", "Genes_count")
for (metric in metrics) {
  log_message("Processing ", metric, "...")
  combined_plot <- create_heatmap(network_stats, metric)
  if (metric == "Edges_count") {
    file_name <- "count_edges"
  } else if (metric == "TFs_count") {
    file_name <- "count_tfs"
  } else if (metric == "Genes_count") {
    file_name <- "count_genes"
  }
  w <- 18
  h <- 10
  ggsave(
    paste0(fig_dir, "/", file_name, ".pdf"),
    combined_plot,
    width = w,
    height = h
  )

  combined_wide <- create_heatmap_wide_region_stage(network_stats, metric)
  n_celltypes <- length(unique(network_stats$Cell_type))

  ggsave(
    paste0(fig_dir, "/", file_name, "_region_stage.pdf"),
    combined_wide,
    width = 20,
    height = 8
  )
}
