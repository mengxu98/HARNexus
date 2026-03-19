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

has_pos_neg <- all(c(
  "Edges_pos_count", "Edges_neg_count",
  "Genes_pos_count", "Genes_neg_count"
) %in% names(network_stats))

if (has_pos_neg) {
  log_message("Creating pos/neg visualizations (scatter, boxplot)...")

  df <- network_stats %>%
    filter(Edges_pos_count + Edges_neg_count > 0) %>%
    mutate(
      pos_frac_edges = Edges_pos_count / (Edges_pos_count + Edges_neg_count),
      pos_frac_genes = ifelse(
        Genes_pos_count + Genes_neg_count > 0,
        Genes_pos_count / (Genes_pos_count + Genes_neg_count),
        NA_real_
      )
    )

  if (nrow(df) > 0) {
    all_regions <- sort(unique(df$Region))
    all_stages <- sort(unique(df$Stage))
    stage_num <- as.numeric(gsub("S", "", all_stages))
    all_stages <- all_stages[order(stage_num)]
    n_regions <- length(all_regions)
    n_stages <- length(all_stages)
    stage_colors <- color_stages[all_stages]
    stage_colors[is.na(stage_colors)] <- "grey50"
    stage_colors <- setNames(stage_colors, all_stages)
    region_shapes <- setNames(
      (seq_len(n_regions) - 1L) %% 25,
      all_regions
    )

    df_plot_edges <- df %>%
      filter(Edges_pos_count > 0, Edges_neg_count > 0)

    pred_edges <- df_plot_edges %>%
      group_by(Cell_type) %>%
      group_modify(~ {
        if (nrow(.x) < 2) {
          return(tibble(Edges_pos_count = numeric(), Edges_neg_count = numeric(), lwr = numeric(), upr = numeric()))
        }
        fit <- lm(
          log10(Edges_neg_count + 1) ~ log10(Edges_pos_count + 1),
          data = .x
        )
        x_r <- range(.x$Edges_pos_count, na.rm = TRUE)
        x_seq <- 10^seq(log10(x_r[1] + 1), log10(x_r[2] + 1), length.out = 50) - 1
        x_seq <- pmax(x_seq, 0.1)
        pred <- predict(
          fit,
          newdata = data.frame(Edges_pos_count = x_seq),
          interval = "confidence"
        )
        tibble(
          Edges_pos_count = x_seq,
          Edges_neg_count = pmax(10^pred[, "fit"] - 1, 0.1),
          lwr = pmax(10^pred[, "lwr"] - 1, 0.1),
          upr = pmax(10^pred[, "upr"] - 1, 0.1)
        )
      }) %>%
      ungroup()

    cor_edges <- df_plot_edges %>%
      group_by(Cell_type) %>%
      summarise(
        r = cor(
          log10(Edges_pos_count + 1),
          log10(Edges_neg_count + 1),
          use = "pairwise.complete.obs"
        ),
        p = tryCatch(
          cor.test(
            log10(Edges_pos_count + 1),
            log10(Edges_neg_count + 1),
            exact = FALSE
          )$p.value,
          error = function(e) NA_real_
        ),
        x_ann = max(min(Edges_pos_count, na.rm = TRUE), 0.1) * 1.05,
        y_ann = max(max(Edges_neg_count, na.rm = TRUE), 1) * 0.95,
        .groups = "drop"
      ) %>%
      mutate(
        label = sprintf(
          "r = %.3f\nP %s",
          r,
          ifelse(is.na(p) | p < 0.001, "< 0.001", sprintf("= %.3f", p))
        )
      )

    cell_types_order <- sort(unique(df$Cell_type))
    strip_labeller <- setNames(
      sprintf(
        "<span style='color:%s'>%s</span>",
        ifelse(
          cell_types_order %in% names(color_celltypes),
          color_celltypes[cell_types_order],
          "black"
        ),
        cell_types_order
      ),
      cell_types_order
    )

    p_scatter_edges <- ggplot(df_plot_edges, aes(
      x = Edges_pos_count,
      y = Edges_neg_count,
      color = factor(Stage, levels = all_stages),
      shape = Region
    )) +
      geom_ribbon(
        data = pred_edges,
        aes(x = Edges_pos_count, ymin = lwr, ymax = upr),
        inherit.aes = FALSE,
        fill = "grey70",
        alpha = 0.3
      ) +
      geom_line(
        data = pred_edges,
        aes(x = Edges_pos_count, y = Edges_neg_count),
        inherit.aes = FALSE,
        color = "grey50",
        linewidth = 1
      ) +
      geom_point(size = 3, alpha = 0.9) +
      geom_text(
        data = cor_edges,
        aes(x = x_ann, y = y_ann, label = label),
        hjust = 0,
        vjust = 1.2,
        size = 3.5,
        inherit.aes = FALSE
      ) +
      facet_wrap(
        ~Cell_type,
        ncol = 3,
        scales = "free",
        labeller = as_labeller(strip_labeller)
      ) +
      scale_x_log10(labels = scales::number, oob = scales::squish) +
      scale_y_log10(labels = scales::number, oob = scales::squish) +
      scale_color_manual(values = stage_colors, na.value = "grey50") +
      scale_shape_manual(values = region_shapes, na.value = 1) +
      theme_bw() +
      theme(
        legend.position = "right",
        strip.background = element_rect(fill = "white"),
        strip.text = ggtext::element_markdown(size = 8),
        aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      labs(
        x = "Positive edge count",
        y = "Negative edge count",
        color = "Stage",
        shape = "Brain region"
      ) +
      guides(
        color = guide_legend(ncol = 5),
        shape = guide_legend(ncol = 2, override.aes = list(size = 3))
      )
    ggsave(
      paste0(fig_dir, "/pos_neg_edges_scatter.pdf"),
      p_scatter_edges,
      width = 12,
      height = 7
    )

    df_genes <- df %>% filter(!is.na(pos_frac_genes))
    df_plot_genes <- df_genes %>%
      filter(Genes_pos_count > 0, Genes_neg_count > 0)
    if (nrow(df_plot_genes) > 0) {
      pred_genes <- df_plot_genes %>%
        group_by(Cell_type) %>%
        group_modify(~ {
          if (nrow(.x) < 2) {
            return(tibble(Genes_pos_count = numeric(), Genes_neg_count = numeric(), lwr = numeric(), upr = numeric()))
          }
          fit <- lm(
            log10(Genes_neg_count + 1) ~ log10(Genes_pos_count + 1),
            data = .x
          )
          x_r <- range(.x$Genes_pos_count, na.rm = TRUE)
          x_seq <- 10^seq(log10(x_r[1] + 1), log10(x_r[2] + 1), length.out = 50) - 1
          x_seq <- pmax(x_seq, 0.1)
          pred <- predict(
            fit,
            newdata = data.frame(Genes_pos_count = x_seq),
            interval = "confidence"
          )
          tibble(
            Genes_pos_count = x_seq,
            Genes_neg_count = pmax(10^pred[, "fit"] - 1, 0.1),
            lwr = pmax(10^pred[, "lwr"] - 1, 0.1),
            upr = pmax(10^pred[, "upr"] - 1, 0.1)
          )
        }) %>%
        ungroup()

      cor_genes <- df_plot_genes %>%
        group_by(Cell_type) %>%
        summarise(
          r = cor(
            log10(Genes_pos_count + 1),
            log10(Genes_neg_count + 1),
            use = "pairwise.complete.obs"
          ),
          p = tryCatch(
            cor.test(
              log10(Genes_pos_count + 1),
              log10(Genes_neg_count + 1),
              exact = FALSE
            )$p.value,
            error = function(e) NA_real_
          ),
          x_ann = max(min(Genes_pos_count, na.rm = TRUE), 0.1) * 1.05,
          y_ann = max(max(Genes_neg_count, na.rm = TRUE), 1) * 0.95,
          .groups = "drop"
        ) %>%
        mutate(
          label = sprintf(
            "r = %.3f\nP %s",
            r,
            ifelse(is.na(p) | p < 0.001, "< 0.001", sprintf("= %.3f", p))
          )
        )

      p_scatter_genes <- ggplot(
        df_plot_genes,
        aes(
          x = Genes_pos_count,
          y = Genes_neg_count,
          color = factor(Stage, levels = all_stages),
          shape = Region
        )
      ) +
        geom_ribbon(
          data = pred_genes,
          aes(x = Genes_pos_count, ymin = lwr, ymax = upr),
          inherit.aes = FALSE,
          fill = "grey70",
          alpha = 0.3
        ) +
        geom_line(
          data = pred_genes,
          aes(x = Genes_pos_count, y = Genes_neg_count),
          inherit.aes = FALSE,
          color = "grey50",
          linewidth = 1
        ) +
        geom_point(size = 3, alpha = 0.9) +
        geom_text(
          data = cor_genes,
          aes(x = x_ann, y = y_ann, label = label),
          hjust = 0,
          vjust = 1.2,
          size = 3.5,
          inherit.aes = FALSE
        ) +
        facet_wrap(
          ~Cell_type,
          ncol = 3,
          scales = "free",
          labeller = as_labeller(strip_labeller)
        ) +
        scale_x_log10(labels = scales::number, oob = scales::squish) +
        scale_y_log10(labels = scales::number, oob = scales::squish) +
        scale_color_manual(values = stage_colors, na.value = "grey50") +
        scale_shape_manual(values = region_shapes, na.value = 1) +
        theme_bw() +
        theme(
          legend.position = "right",
          strip.background = element_rect(fill = "white"),
          strip.text = ggtext::element_markdown(size = 8),
          aspect.ratio = 1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        ) +
        labs(
          x = "Positively HAR target gene count",
          y = "Negatively HAR target gene count",
          color = "Stage",
          shape = "Brain region",
        ) +
        guides(
          color = guide_legend(ncol = 5),
          shape = guide_legend(ncol = 2, override.aes = list(size = 3))
        )
      ggsave(
        paste0(fig_dir, "/pos_neg_genes_scatter.pdf"),
        p_scatter_genes,
        width = 11.5,
        height = 7
      )

      p_scatter <- p_scatter_edges + p_scatter_genes + plot_layout(
        ncol = 2, guides = "collect"
      )
      ggsave(
        paste0(fig_dir, "/pos_neg_scatter.pdf"),
        p_scatter,
        width = 20,
        height = 7
      )
    }

    cell_types_box <- sort(unique(df$Cell_type))
    mean_df <- df %>%
      group_by(Cell_type) %>%
      summarise(
        mean_val = mean(pos_frac_edges, na.rm = TRUE), .groups = "drop"
      ) %>%
      mutate(
        pos = match(Cell_type, cell_types_box),
        xmin = pos - 0.45,
        xmax = pos + 0.45,
        fill_leg = "Mean"
      )

    p_box <- ggplot(
      df,
      aes(x = Cell_type, y = pos_frac_edges, fill = Cell_type)
    ) +
      # geom_rect(
      #   data = mean_df,
      #   aes(
      #     xmin = xmin,
      #     xmax = xmax,
      #     ymin = 0,
      #     ymax = 1
      #   ),
      #   fill = "grey90",
      #   inherit.aes = FALSE
      # ) +
      geom_rect(
        data = mean_df,
        aes(
          xmin = xmin,
          xmax = xmax,
          ymin = 0,
          ymax = mean_val,
          fill = fill_leg
        ),
        # alpha = 0.5,
        inherit.aes = FALSE
      ) +
      geom_boxplot(outlier.alpha = 0.5, width = 0.7, show.legend = FALSE) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey10") +
      scale_x_discrete(limits = cell_types_box) +
      scale_fill_manual(
        values = c(color_celltypes, "Mean" = "grey60"),
        guide = guide_legend(
          override.aes = list(alpha = 1),
          title = NULL
        ),
        breaks = "Mean"
      ) +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      labs(
        x = "Cell type",
        y = "Proportion of positive edges"
      )
    ggsave(
      paste0(fig_dir, "/pos_neg_edges_proportion_boxplot.pdf"),
      p_box,
      width = 4,
      height = 3.5
    )
  }
} else {
  log_message(
    "Skipping pos/neg visualizations: run 04_altas_statistics.py first"
  )
}
