plot_auc_dimred <- function(
    seurat,
    gene_set_name,
    reduction = "umap",
    point_size = 0.5,
    highlight_celltypes = NULL,
    highlight_color = "#D70440",
    show_auc_legend = FALSE,
    color_methods) {
  dim_coords <- as.data.frame(
    seurat@reductions[[reduction]]@cell.embeddings
  )

  plot_data <- data.frame(
    Dim1 = dim_coords[, 1],
    Dim2 = dim_coords[, 2],
    AUC = as.numeric(unlist(seurat[[paste0("AUC_", gene_set_name)]])),
    Pass_threshold = as.logical(
      unlist(
        seurat[[paste0("Pass_threshold_", gene_set_name)]]
      )
    ),
    CellType = seurat$CellType
  )

  centroids <- plot_data %>%
    group_by(CellType) %>%
    summarise(
      Centroid_x = mean(Dim1),
      Centroid_y = mean(Dim2)
    )

  dim_labels <- c("UMAP_1", "UMAP_2")
  colnames(plot_data)[1:2] <- dim_labels

  method_name <- strsplit(gene_set_name, " ")[[1]][1]
  if (grepl("Random", method_name)) {
    method_name <- gene_set_name
  }
  gene_set_color <- color_methods[method_name]

  n_breaks <- 1000
  color_pal_neg <- colorRampPalette(
    c("black", "blue", "skyblue")
  )(n_breaks / 2)
  color_pal_pos <- colorRampPalette(
    c("pink", "magenta", "red")
  )(n_breaks / 2)
  full_palette <- c(color_pal_neg, color_pal_pos)

  p <- ggplot(
    plot_data,
    aes(x = .data[[dim_labels[1]]], y = .data[[dim_labels[2]]])
  ) +
    ggrastr::geom_point_rast(
      aes(color = .data$AUC),
      size = 0.5, raster.dpi = 300
    ) +
    scale_color_gradientn(
      colors = full_palette,
      limits = range(plot_data$AUC, na.rm = TRUE),
      guide = if (show_auc_legend) "colorbar" else "none"
    ) +
    labs(
      title = gene_set_name,
      x = dim_labels[1],
      y = dim_labels[2],
      color = if (show_auc_legend) "AUC" else NULL
    ) +
    coord_fixed(ratio = 1) +
    theme_classic() +
    theme(
      plot.title = element_text(
        color = gene_set_color
      ),
      legend.position = if (show_auc_legend) "bottom" else "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(1, 1, 1, 1, "pt")
    )

  if (!is.null(highlight_celltypes) && length(highlight_celltypes) > 0) {
    for (celltype in highlight_celltypes) {
      celltype_data <- plot_data[plot_data$CellType == celltype, ]
      if (nrow(celltype_data) > 0) {
        p <- p + stat_ellipse(
          data = celltype_data,
          aes(x = .data[[dim_labels[1]]], y = .data[[dim_labels[2]]]),
          type = "t",
          level = 0.98,
          color = highlight_color,
          linewidth = 0.6,
          linetype = "3313",
          inherit.aes = FALSE
        )

        celltype_centroid <- centroids[centroids$CellType == celltype, ]
        if (nrow(celltype_centroid) > 0) {
          x_data <- celltype_data[[dim_labels[1]]]
          y_data <- celltype_data[[dim_labels[2]]]

          x_center <- mean(x_data, na.rm = TRUE)
          y_center <- mean(y_data, na.rm = TRUE)

          cov_matrix <- cov(cbind(x_data, y_data), use = "complete.obs")
          x_sd <- sqrt(cov_matrix[1, 1])
          n <- length(x_data)
          df <- n - 1
          f_crit <- qf(0.98, 2, df)
          scale_factor <- sqrt(f_crit * 2)
          label_x <- x_center + scale_factor * x_sd
          label_y <- y_center

          p <- p + annotate(
            "text",
            x = label_x,
            y = label_y,
            label = celltype,
            color = highlight_color,
            size = 3.5,
            hjust = 0,
            vjust = 0.5
          )
        }
      }
    }
  }

  return(p)
}

prepare_auc_data <- function(
    seurat,
    selected_thresholds,
    color_methods) {
  auc_data_list <- lapply(
    names(selected_thresholds),
    function(gene_set_name) {
      method_name <- strsplit(gene_set_name, " ")[[1]][1]
      if (grepl("Random", method_name)) {
        method_name <- gene_set_name
      }

      data.frame(
        GeneSet = gene_set_name,
        Method = method_name,
        AUC = as.numeric(unlist(seurat[[paste0("AUC_", gene_set_name)]])),
        CellType = seurat$CellType
      )
    }
  )
  result <- do.call(rbind, auc_data_list)

  result$Method <- factor(
    result$Method,
    levels = names(color_methods)
  )
  result$GeneSet <- factor(
    result$GeneSet,
    levels = names(selected_thresholds)
  )

  return(result)
}

threshold_proportion_data <- function(
    seurat,
    selected_thresholds,
    color_methods) {
  threshold_data <- lapply(
    names(selected_thresholds),
    function(gene_set_name) {
      method_name <- strsplit(gene_set_name, " ")[[1]][1]
      if (grepl("Random", method_name)) {
        method_name <- gene_set_name
      }

      col_name <- paste0("Pass_threshold_", gene_set_name)
      pass_threshold <- as.logical(unlist(seurat[[col_name]]))
      cell_types <- seurat$CellType

      result <- data.frame(
        GeneSet = gene_set_name,
        Method = method_name,
        CellType = factor(cell_types),
        Status = factor(
          ifelse(
            pass_threshold, "Above AUC threshold", "Below AUC threshold"
          )
        )
      )

      return(result)
    }
  )

  threshold_data <- do.call(rbind, threshold_data)

  threshold_data <- threshold_data %>%
    group_by(Method, GeneSet, CellType, Status) %>%
    summarise(Count = n(), .groups = "drop") %>%
    group_by(Method, GeneSet, CellType) %>%
    mutate(Proportion = Count / sum(Count))

  threshold_data$Method <- factor(
    threshold_data$Method,
    levels = names(color_methods)
  )

  return(threshold_data)
}

aucell_analysis <- function(
    seurat,
    gene_sets_list,
    random_sets = TRUE,
    random_sizes = c(100, 300),
    output_dir = "",
    seed = 2025,
    reduction = "umap",
    threshold_offset = 0,
    use_adaptive_threshold = FALSE,
    adaptive_quantile_base = 0.85,
    adaptive_quantile_range = 0.15,
    point_size = 0.5,
    boxplot_y_range = TRUE,
    color_celltypes = c(
      "Radial glia" = "#8076A3",
      "Neuroblasts" = "#ED5736",
      "Excitatory neurons" = "#0AA344",
      "Inhibitory neurons" = "#2177B8",
      "Astrocytes" = "#D70440",
      "Oligodendrocyte progenitor cells" = "#F9BD10",
      "Oligodendrocytes" = "#B14B28",
      "Microglia" = "#006D87",
      "Endothelial cells" = "#5E7987"
    ),
    color_methods = c(
      "GENIE3" = "#2177B8",
      "HARNexus" = "#D70440",
      "LEAP" = "#F9BD10",
      "PPCOR" = "#0AA344",
      "GREAT" = "#8076A3",
      "TFs" = "#006D87",
      "All" = "#ED5736",
      "Random (100)" = "#ee8484",
      "Random (300)" = "#d83e2a"
    ),
    highlight_celltypes = NULL,
    highlight_color = "#D70440",
    width = 13,
    height = 12) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  set.seed(seed)

  cell_counts <- table(seurat$CellType)
  seurat$cluster_num <- paste0(
    seurat$CellType,
    "\n(",
    cell_counts[match(seurat$CellType, names(cell_counts))],
    ")"
  )

  expr_matrix <- as.matrix(
    GetAssayData(
      seurat,
      layer = "data"
    )
  )

  gene_sets_collection <- list()
  for (name in names(gene_sets_list)) {
    gene_set <- gene_sets_list[[name]]
    gene_set <- gene_set[gene_set %in% rownames(expr_matrix)]
    gene_sets_collection[[name]] <- GSEABase::GeneSet(
      gene_set,
      setName = paste0(name, " (", length(gene_set), ")")
    )
  }

  if (random_sets) {
    all_genes <- NULL
    if ("All" %in% names(gene_sets_list)) {
      all_genes <- gene_sets_list[["All"]]
      all_genes <- all_genes[all_genes %in% rownames(expr_matrix)]
    }

    if (!is.null(all_genes) && length(all_genes) > 0) {
      for (size in random_sizes) {
        actual_size <- min(size, length(all_genes))
        gene_sets_collection[[paste0("Random (", size, ")")]] <- GSEABase::GeneSet(
          sample(all_genes, actual_size),
          setName = paste0("Random (", size, ")")
        )
      }
    }
  }

  gene_sets <- GSEABase::GeneSetCollection(
    gene_sets_collection
  )

  cells_rankings <- AUCell::AUCell_buildRankings(
    expr_matrix,
    plotStats = FALSE
  )
  cells_auc <- AUCell::AUCell_calcAUC(
    gene_sets,
    cells_rankings
  )

  pdf(
    file = file.path(output_dir, "hist_plots.pdf"),
    width = 8,
    height = 3.5
  )
  par(
    mfrow = c(2, 4),
    mar = c(2, 2, 2, 1),
    cex.main = 1,
    cex.lab = 1,
    cex.axis = 1
  )
  cells_assignment <- AUCell::AUCell_exploreThresholds(
    cells_auc,
    plotHist = TRUE,
    assignCells = TRUE
  )
  dev.off()

  selected_thresholds <- AUCell::getThresholdSelected(cells_assignment)

  auc_data <- AUCell::getAUC(cells_auc)

  auc_stats <- lapply(names(selected_thresholds), function(geneSetName) {
    auc_values <- as.numeric(auc_data[geneSetName, ])
    list(
      mean = mean(auc_values, na.rm = TRUE),
      median = median(auc_values, na.rm = TRUE),
      q90 = quantile(auc_values, 0.9, na.rm = TRUE)
    )
  })
  names(auc_stats) <- names(selected_thresholds)

  adaptive_quantiles <- NULL
  if (use_adaptive_threshold) {
    all_means <- sapply(auc_stats, function(x) x$mean)
    mean_min <- min(all_means, na.rm = TRUE)
    mean_max <- max(all_means, na.rm = TRUE)
    mean_range <- mean_max - mean_min

    adaptive_quantiles <- numeric(length(selected_thresholds))
    names(adaptive_quantiles) <- names(selected_thresholds)

    for (geneSetName in names(selected_thresholds)) {
      current_mean <- auc_stats[[geneSetName]]$mean

      if (mean_range > 0) {
        normalized_mean <- (current_mean - mean_min) / mean_range
      } else {
        normalized_mean <- 0.5
      }
      quantile_adjustment <- normalized_mean * adaptive_quantile_range
      adaptive_quantile <- adaptive_quantile_base - quantile_adjustment

      adaptive_quantile <- max(0.5, min(0.95, adaptive_quantile))
      adaptive_quantiles[geneSetName] <- adaptive_quantile

      auc_values <- as.numeric(auc_data[geneSetName, ])
      selected_thresholds[geneSetName] <- quantile(
        auc_values,
        probs = adaptive_quantile,
        na.rm = TRUE
      )
    }
  }

  threshold_summary <- data.frame(
    GeneSet = character(),
    ThresholdType = character(),
    Threshold = numeric(),
    AdaptiveQuantile = numeric(),
    AUC_Mean = numeric(),
    AUC_Median = numeric(),
    AUC_Min = numeric(),
    AUC_Max = numeric(),
    AUC_Q25 = numeric(),
    AUC_Q75 = numeric(),
    AUC_Q90 = numeric(),
    Cells_Above = integer(),
    Cells_Total = integer(),
    Proportion_Above = numeric(),
    stringsAsFactors = FALSE
  )

  for (geneSetName in names(selected_thresholds)) {
    auc_values <- as.numeric(auc_data[geneSetName, ])
    base_threshold <- selected_thresholds[geneSetName]
    adjusted_threshold <- base_threshold - threshold_offset

    pass_threshold <- auc_values > adjusted_threshold
    n_above <- sum(pass_threshold, na.rm = TRUE)
    n_total <- length(auc_values)

    threshold_type <- if (use_adaptive_threshold) "Adaptive threshold" else "Automatic threshold"
    adaptive_q <- if (use_adaptive_threshold && !is.null(adaptive_quantiles)) {
      adaptive_quantiles[geneSetName]
    } else {
      NA
    }

    threshold_summary <- rbind(
      threshold_summary,
      data.frame(
        GeneSet = geneSetName,
        ThresholdType = threshold_type,
        Threshold = base_threshold,
        AdaptiveQuantile = adaptive_q,
        AUC_Mean = mean(auc_values, na.rm = TRUE),
        AUC_Median = median(auc_values, na.rm = TRUE),
        AUC_Min = min(auc_values, na.rm = TRUE),
        AUC_Max = max(auc_values, na.rm = TRUE),
        AUC_Q25 = quantile(auc_values, 0.25, na.rm = TRUE),
        AUC_Q75 = quantile(auc_values, 0.75, na.rm = TRUE),
        AUC_Q90 = quantile(auc_values, 0.90, na.rm = TRUE),
        Cells_Above = n_above,
        Cells_Total = n_total,
        Proportion_Above = n_above / n_total,
        stringsAsFactors = FALSE
      )
    )
  }

  for (geneSetName in names(selected_thresholds)) {
    auc_values <- as.numeric(auc_data[geneSetName, ])
    base_threshold <- selected_thresholds[geneSetName]
    adjusted_threshold <- base_threshold - threshold_offset
    pass_threshold <- auc_values > adjusted_threshold

    cell_types <- seurat$CellType
    cell_type_stats <- table(cell_types[pass_threshold])
    cell_type_total <- table(cell_types)

    for (ct in names(cell_type_total)) {
      n_above <- ifelse(ct %in% names(cell_type_stats), cell_type_stats[ct], 0)
      n_total <- cell_type_total[ct]
      prop <- n_above / n_total * 100
    }
  }
  if (output_dir != "") {
    write.csv(
      threshold_summary,
      file = file.path(output_dir, "threshold_summary.csv"),
      row.names = FALSE
    )
  }

  for (geneSetName in names(selected_thresholds)) {
    seurat[[paste0("AUC_", geneSetName)]] <- auc_data[geneSetName, ]
    seurat[[paste0("Pass_threshold_", geneSetName)]] <-
      auc_data[geneSetName, ] > (selected_thresholds[geneSetName] - threshold_offset)
  }

  auc_stats_data <- prepare_auc_data(
    seurat,
    selected_thresholds,
    color_methods
  )
  threshold_data <- threshold_proportion_data(
    seurat,
    selected_thresholds,
    color_methods
  )

  celltype_levels <- sort(as.character(unique(seurat$CellType)))
  auc_stats_data$CellType <- factor(
    as.character(auc_stats_data$CellType),
    levels = celltype_levels
  )
  threshold_data$CellType <- factor(
    as.character(threshold_data$CellType),
    levels = celltype_levels
  )

  if (boxplot_y_range) {
    y_axis_range <- range(auc_stats_data$AUC, na.rm = TRUE)
    y_axis_range[1] <- max(0, y_axis_range[1] - 0.05)
    y_axis_range[2] <- y_axis_range[2] + 0.05
  } else {
    y_axis_range <- NULL
  }

  method_plots <- lapply(
    seq_along(levels(auc_stats_data$GeneSet)),
    function(plot_idx) {
      gene_set_name <- levels(auc_stats_data$GeneSet)[plot_idx]
      is_top_row <- plot_idx <= 4
      method_name <- strsplit(gene_set_name, " ")[[1]][1]
      if (grepl("Random", method_name)) {
        method_name <- gene_set_name
      }

      subset_data <- auc_stats_data[auc_stats_data$GeneSet == gene_set_name, ]
      subset_data <- subset_data %>%
        group_by(CellType) %>%
        sample_frac(0.1) %>%
        ungroup()

      p_box <- ggplot(subset_data, aes(x = CellType, y = AUC)) +
        geom_boxplot(
          aes(fill = CellType),
          width = 0.9,
          outlier.shape = NA
        ) +
        geom_jitter(
          width = 0.3,
          size = 0.2,
          alpha = 0.2
        ) +
        scale_fill_manual(
          values = color_celltypes,
          guide = "none"
        ) +
        theme_bw() +
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          plot.margin = margin(-5, 5, -5, 5)
        ) +
        labs(
          x = NULL,
          y = "AUC"
        )

      current_y_range <- if (boxplot_y_range && !is.null(y_axis_range)) {
        y_axis_range
      } else {
        range(subset_data$AUC, na.rm = TRUE)
      }

      text_y <- NULL
      if (!is.null(highlight_celltypes) && length(highlight_celltypes) > 0) {
        highlight_data <- subset_data[subset_data$CellType %in% highlight_celltypes, ]
        if (nrow(highlight_data) > 0) {
          mean_auc <- mean(highlight_data$AUC, na.rm = TRUE)
          celltype_levels <- levels(factor(subset_data$CellType))
          highlight_positions <- which(celltype_levels %in% highlight_celltypes)

          text_y <- current_y_range[2] * 0.9
          p_box <- p_box + geom_hline(
            yintercept = mean_auc,
            color = highlight_color,
            linetype = "dashed"
          )

          for (pos in highlight_positions) {
            p_box <- p_box + annotate(
              "text",
              x = pos,
              y = text_y,
              label = sprintf("Mean AUC: %.3f", mean_auc),
              color = highlight_color,
              size = 3,
              hjust = 0.5,
              vjust = 0
            )
          }
        }
      }

      if (boxplot_y_range) {
        p_box <- p_box + scale_y_continuous(
          limits = current_y_range,
          expand = c(0, 0),
          labels = function(x) ifelse(x == 0, "", x)
        )
      } else {
        p_box <- p_box + scale_y_continuous(
          expand = c(0, 0.05),
          labels = function(x) ifelse(x == 0, "", x)
        )
      }

      method_data <- threshold_data[threshold_data$GeneSet == gene_set_name, ]
      method_data$Status <- factor(
        method_data$Status,
        levels = c("Above AUC threshold", "Below AUC threshold")
      )

      summary_data <- data.frame(
        CellType = unique(method_data$CellType)
      )

      if (nrow(method_data) > 0) {
        total_counts <- aggregate(
          Count ~ CellType,
          data = method_data,
          sum
        )
        summary_data$Total <- total_counts$Count[match(
          summary_data$CellType,
          total_counts$CellType
        )]

        above_data <- method_data[method_data$Status == "Above AUC threshold", ]
        if (nrow(above_data) > 0) {
          above_counts <- aggregate(
            Count ~ CellType,
            data = above_data,
            sum
          )
          summary_data$Above <- above_counts$Count[match(
            summary_data$CellType,
            above_counts$CellType
          )]
        } else {
          summary_data$Above <- 0
        }
      } else {
        summary_data$Total <- 0
        summary_data$Above <- 0
      }

      summary_data$Total[is.na(summary_data$Total)] <- 0
      summary_data$Above[is.na(summary_data$Above)] <- 0

      summary_data$Percent <- round(
        summary_data$Above / summary_data$Total * 100, 1
      )
      summary_data$label <- sprintf(
        "%d/%d\n(%.1f%%)",
        summary_data$Above,
        summary_data$Total,
        summary_data$Percent
      )

      p_threshold <- ggplot() +
        geom_bar(
          data = method_data,
          aes(x = CellType, y = Proportion, fill = Status),
          stat = "identity",
          position = "stack",
          width = 0.9
        ) +
        scale_fill_manual(
          values = c(
            "Above AUC threshold" = "black",
            "Below AUC threshold" = "gray50"
          ),
          name = "Status"
        ) +
        geom_text(
          data = summary_data,
          aes(x = CellType, y = 0.5, label = label),
          angle = 90,
          size = 2.8,
          color = "white",
          vjust = 0.5,
          hjust = 0.5
        ) +
        theme_bw() +
        theme(
          axis.text.x = if (is_top_row) element_blank() else element_text(angle = 30, hjust = 1),
          axis.ticks.x = if (is_top_row) element_blank() else element_line(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          plot.margin = margin(-5, 5, 10, 5)
        ) +
        labs(x = if (is_top_row) NULL else "Celltype", y = "Proportion")

      combined_plot <- p_box / p_threshold +
        plot_layout(heights = c(0.6, 0.4)) &
        theme(panel.spacing = unit(0, "mm"))

      return(combined_plot)
    }
  )

  auc_plots <- lapply(
    names(selected_thresholds), function(gene_set_name) {
      method_name <- strsplit(gene_set_name, " ")[[1]][1]
      if (grepl("Random", method_name)) {
        method_name <- gene_set_name
      }
      show_auc_legend <- (method_name == "HARNexus")

      plot_auc_dimred(
        seurat,
        gene_set_name,
        reduction = reduction,
        point_size = 0.5,
        color_methods = color_methods,
        highlight_celltypes = highlight_celltypes,
        highlight_color = highlight_color,
        show_auc_legend = show_auc_legend
      )
    }
  )
  names(auc_plots) <- names(selected_thresholds)

  dim_coords <- as.data.frame(
    seurat@reductions[[reduction]]@cell.embeddings
  )
  dim_labels <- c("UMAP_1", "UMAP_2")

  celltype_plot_data <- data.frame(
    Dim1 = dim_coords[, 1],
    Dim2 = dim_coords[, 2],
    CellType = seurat$CellType
  )
  colnames(celltype_plot_data)[1:2] <- dim_labels

  celltype_umap_plot <- ggplot(
    celltype_plot_data,
    aes(x = .data[[dim_labels[1]]], y = .data[[dim_labels[2]]])
  ) +
    ggrastr::geom_point_rast(
      aes(color = .data$CellType),
      size = 0.5, raster.dpi = 300
    ) +
    scale_color_manual(values = color_celltypes, guide = "none") +
    labs(
      title = "Celltype",
      x = dim_labels[1],
      y = dim_labels[2],
      color = "Celltype"
    ) +
    coord_fixed(ratio = 1) +
    theme_classic() +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(1, 1, 1, 1, "pt")
    )

  actual_celltypes <- unique(celltype_plot_data$CellType)
  actual_celltypes <- actual_celltypes[!is.na(actual_celltypes)]
  legend_colors <- color_celltypes[names(color_celltypes) %in% actual_celltypes]
  legend_colors <- legend_colors[order(match(names(legend_colors), actual_celltypes))]

  legend_labels <- names(legend_colors)
  legend_labels[legend_labels == "Oligodendrocyte progenitor cells"] <- "Oligodendrocyte\nprogenitor cells"

  legend_data <- data.frame(
    CellType = factor(names(legend_colors), levels = names(legend_colors)),
    CellTypeLabel = legend_labels,
    x = 0.3,
    y = seq_along(legend_colors)
  )

  legend_plot <- ggplot(legend_data, aes(x = x, y = y)) +
    geom_point(
      aes(color = CellType),
      size = 4,
      shape = 15
    ) +
    scale_color_manual(values = legend_colors, guide = "none") +
    geom_text(
      aes(label = CellTypeLabel),
      x = 0.5,
      hjust = 0,
      vjust = 0.5,
      size = 3.2
    ) +
    annotate(
      "text",
      x = 0.4,
      y = length(legend_colors) + 0.8,
      label = "Celltype",
      hjust = 0,
      vjust = 0,
      size = 4
    ) +
    xlim(0.2, 2.5) +
    ylim(0.5, length(legend_colors) + 1.2) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin = margin(-5, 5, 10, 5)
    )

  # layout of 5 columns and 4 rows:
  # first row:  [Cell type UMAP] + [UMAP1] + [UMAP2] + [UMAP3] + [UMAP4]
  # second row: [Legend] +         [Stat1] + [Stat2] + [Stat3] + [Stat4]
  # third row:  [UMAP5] +          [UMAP6] + [UMAP7] + [UMAP8] + [UMAP9]
  # fourth row: [Stat5] +          [Stat6] + [Stat7] + [Stat8] + [Stat9]

  final_plot_all <- celltype_umap_plot +
    auc_plots[[1]] +
    auc_plots[[2]] +
    auc_plots[[3]] +
    auc_plots[[4]] +
    legend_plot +
    method_plots[[1]] +
    method_plots[[2]] +
    method_plots[[3]] +
    method_plots[[4]] +
    auc_plots[[5]] +
    auc_plots[[6]] +
    auc_plots[[7]] +
    auc_plots[[8]] +
    auc_plots[[9]] +
    method_plots[[5]] +
    method_plots[[6]] +
    method_plots[[7]] +
    method_plots[[8]] +
    method_plots[[9]] +
    plot_layout(
      guides = "collect",
      ncol = 5,
      nrow = 4
    ) & theme(
    legend.position = "bottom",
    plot.margin = margin(1, 1, 1, 1, "pt")
  ) &
    guides(shape = "none", linetype = "none")

  ggsave(
    file.path(output_dir, "stats_plots.pdf"),
    final_plot_all,
    width = width,
    height = height
  )
  return(seurat)
}
