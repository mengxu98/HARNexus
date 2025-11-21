source("code/functions/utils.R")


plot_auc_dimred <- function(
    seurat,
    gene_set_name,
    reduction = "umap",
    point_size = 0.5,
    show_labels = TRUE,
    color_palette_method) {
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
    CellType = seurat$cell_name
  )

  centroids <- plot_data %>%
    group_by(CellType) %>%
    summarise(
      Centroid_x = mean(Dim1),
      Centroid_y = mean(Dim2)
    )

  dim_labels <- if (reduction == "umap") {
    c("UMAP_1", "UMAP_2")
  } else {
    c("t-SNE_1", "t-SNE_2")
  }
  colnames(plot_data)[1:2] <- dim_labels

  method_name <- strsplit(gene_set_name, " ")[[1]][1]
  if (grepl("Random", method_name)) {
    method_name <- gene_set_name
  }
  gene_set_color <- color_palette_method[method_name]

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
    geom_point(aes(color = .data$AUC), size = point_size) +
    scale_color_gradientn(
      colors = full_palette,
      limits = range(plot_data$AUC, na.rm = TRUE)
    ) +
    labs(
      title = gene_set_name,
      x = dim_labels[1],
      y = dim_labels[2],
      color = "AUC"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        color = gene_set_color,
        face = "bold"
      ),
      legend.position = "none",
      plot.margin = margin(5, 5, -5, 5)
    )

  if (show_labels) {
    p <- p + geom_label(
      data = centroids,
      aes(x = Centroid_x, y = Centroid_y, label = CellType),
      size = 3.5,
      fill = "white",
      alpha = 0.85,
      label.padding = unit(0.2, "lines"),
      label.size = 0,
      check_overlap = TRUE
    )
  }

  return(p)
}

prepare_auc_data <- function(
    seurat,
    selected_thresholds,
    color_palette_method) {
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
        CellType = seurat$cell_name
      )
    }
  )
  result <- do.call(rbind, auc_data_list)

  result$Method <- factor(
    result$Method,
    levels = names(color_palette_method)
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
    color_palette_method) {
  threshold_data <- lapply(
    names(selected_thresholds),
    function(gene_set_name) {
      method_name <- strsplit(gene_set_name, " ")[[1]][1]
      if (grepl("Random", method_name)) {
        method_name <- gene_set_name
      }

      col_name <- paste0("Pass_threshold_", gene_set_name)
      pass_threshold <- as.logical(unlist(seurat[[col_name]]))
      cell_types <- seurat$cell_name

      result <- data.frame(
        GeneSet = gene_set_name,
        Method = method_name,
        CellType = factor(cell_types),
        Status = factor(
          ifelse(
            pass_threshold, "Above threshold", "Below threshold"
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
    levels = names(color_palette_method)
  )

  return(threshold_data)
}

reduction_function <- function(
    reduction_data,
    selected_thresholds,
    cells_auc,
    threshold_offset) {
  par(
    mfrow = c(2, 4),
    mar = c(2, 2, 2, 1),
    cex.main = 1,
    cex.lab = 1, cex.axis = 1
  )
  for (geneSetName in names(selected_thresholds)) {
    n_breaks <- 1000
    color_pal_neg <- colorRampPalette(
      c("black", "blue", "skyblue")
    )(n_breaks)
    color_pal_pos <- colorRampPalette(
      c("pink", "magenta", "red")
    )(n_breaks)

    pass_threshold <- getAUC(
      cells_auc
    )[geneSetName, ] > (selected_thresholds[geneSetName] - threshold_offset)
    if (sum(pass_threshold) > 0) {
      auc_split <- split(getAUC(cells_auc)[geneSetName, ], pass_threshold)

      cell_color <- c(
        setNames(
          color_pal_neg[cut(auc_split[[1]], breaks = n_breaks)],
          names(auc_split[[1]])
        ),
        setNames(
          color_pal_pos[cut(auc_split[[2]], breaks = n_breaks)],
          names(auc_split[[2]])
        )
      )

      plot(
        reduction_data,
        main = geneSetName,
        sub = "Pink/red cells pass the threshold",
        col = cell_color[rownames(reduction_data)],
        pch = 16
      )
    }
  }
}

aucell_analysis <- function(
    seurat,
    gene_sets_list,
    var_genes = NULL,
    random_sets = TRUE,
    random_sizes = c(100, 200, 400),
    output_dir = "",
    seed = 2024,
    reductions = c("tsne", "umap"),
    threshold_offset = 0.02,
    point_size = 0.5,
    color_palette_cluster = c(
      "Astro" = "#005ea3",
      "Endo" = "#24B700",
      "Micro" = "#00C1AB",
      "OPC" = "#00ACFC",
      "ExN" = "#f0749d",
      "InN" = "#c21b90",
      "NPC" = "#e29828",
      "Olig" = "#5865d3",
      "Perc" = "#c08f09"
    ),
    color_palette_method = c(
      "GENIE3" = "#cc3366",
      "LEAP" = "#ff9900",
      "PPCOR" = "#339999",
      "HARNexus" = "#005ea3",
      "GREAT" = "#00a323",
      "Random (100)" = "#ee8484",
      "Random (200)" = "#e34949",
      "Random (400)" = "#d82a2a"
    )) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  set.seed(seed)

  cell_counts <- table(seurat$cell_name)
  seurat$cluster_num <- paste0(
    seurat$cell_name,
    "\n(",
    cell_counts[match(seurat$cell_name, names(cell_counts))],
    ")"
  )
  color_palette_cluster_seurat <- setNames(
    color_palette_cluster,
    paste0(
      names(color_palette_cluster),
      "\n(",
      cell_counts[names(color_palette_cluster)],
      ")"
    )
  )
  tsne_data <- seurat@reductions$tsne@cell.embeddings
  umap_data <- seurat@reductions$umap@cell.embeddings

  expr_matrix <- as.matrix(
    GetAssayData(
      seurat,
      layer = "data"
    )
  )

  if (is.null(var_genes)) {
    var_genes <- VariableFeatures(seurat)
  }

  gene_sets_collection <- list()
  for (name in names(gene_sets_list)) {
    gene_set <- gene_sets_list[[name]]
    gene_set <- gene_set[gene_set %in% rownames(expr_matrix)]
    gene_sets_collection[[name]] <- GSEABase::GeneSet(
      gene_set,
      setName = paste0(name, " (", length(gene_set), ")")
    )
  }

  if (random_sets && !is.null(var_genes)) {
    for (size in random_sizes) {
      gene_sets_collection[[paste0("Random (", size, ")")]] <- GSEABase::GeneSet(
        sample(var_genes, size),
        setName = paste0("Random (", size, ")")
      )
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
    file = paste0(output_dir, "/AUCell_hist_plots.pdf"),
    width = 6.5,
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

  png(
    file = paste0(output_dir, "/AUCell_hist_plots.png"),
    width = 4000,
    height = 4000 * 3.5 / 6.5,
    res = 600
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

  if ("tsne" %in% reductions) {
    pdf(
      file = file.path(output_dir, "AUCell_tsne.pdf"),
      width = 6.5,
      height = 3.5
    )
    reduction_function(
      tsne_data,
      selected_thresholds,
      cells_auc,
      threshold_offset
    )
    dev.off()

    png(
      file = file.path(output_dir, "AUCell_tsne.png"),
      width = 4000,
      height = 4000 * 3.5 / 6.5,
      res = 600
    )
    reduction_function(
      tsne_data,
      selected_thresholds,
      cells_auc,
      threshold_offset
    )
    dev.off()

    pdf(
      file = file.path(output_dir, "AUCell_hist_tsne.pdf"),
      width = 8,
      height = 6
    )
    par(
      mfrow = c(4, 6),
      mar = c(2, 2, 2, 1),
      cex.main = 1,
      cex.lab = 1,
      cex.axis = 1
    )
    AUCell_plotTSNE(
      tSNE = tsne_data,
      exprMat = expr_matrix,
      cellsAUC = cells_auc,
      thresholds = selected_thresholds,
      reorderGeneSets = FALSE,
      cex = 1,
      alphaOn = 1,
      alphaOff = 0.2,
      borderColor = adjustcolor("black", alpha.f = 0),
      offColor = "gray80",
      plots = c("histogram", "binaryAUC", "AUC", "expression"),
      exprCols = c("goldenrod1", "darkorange", "brown"),
      asPNG = FALSE
    )
    dev.off()
  }

  if ("umap" %in% reductions) {
    pdf(
      file = file.path(output_dir, "AUCell_umap.pdf"),
      width = 6.5,
      height = 3.5
    )
    reduction_function(
      umap_data,
      selected_thresholds,
      cells_auc,
      threshold_offset
    )
    dev.off()

    png(
      file = file.path(output_dir, "AUCell_umap.png"),
      width = 4000,
      height = 4000 * 3.5 / 6.5,
      res = 600
    )
    reduction_function(
      umap_data,
      selected_thresholds,
      cells_auc,
      threshold_offset
    )
    dev.off()

    pdf(
      file = file.path(output_dir, "AUCell_hist_umap.pdf"),
      width = 8,
      height = 6
    )
    par(
      mfrow = c(4, 6),
      mar = c(2, 2, 2, 1),
      cex.main = 1,
      cex.lab = 1,
      cex.axis = 1
    )
    AUCell::AUCell_plotTSNE(
      tSNE = umap_data,
      exprMat = expr_matrix,
      cellsAUC = cells_auc,
      thresholds = selected_thresholds,
      reorderGeneSets = FALSE,
      cex = 1,
      alphaOn = 1,
      alphaOff = 0.2,
      borderColor = adjustcolor("black", alpha.f = 0),
      offColor = "gray80",
      plots = c("histogram", "binaryAUC", "AUC", "expression"),
      exprCols = c("goldenrod1", "darkorange", "brown"),
      asPNG = FALSE
    )
    dev.off()
  }

  auc_data <- getAUC(cells_auc)
  for (geneSetName in names(selected_thresholds)) {
    seurat[[paste0("AUC_", geneSetName)]] <- auc_data[geneSetName, ]
    seurat[[paste0("Pass_threshold_", geneSetName)]] <-
      auc_data[geneSetName, ] > (selected_thresholds[geneSetName] - threshold_offset)
  }

  for (reduction in reductions) {
    auc_plots <- lapply(
      names(selected_thresholds), function(gene_set_name) {
        plot_auc_dimred(
          seurat,
          gene_set_name,
          reduction = reduction,
          show_labels = FALSE,
          point_size = 0.5,
          color_palette_method
        )
      }
    )
    names(auc_plots) <- names(selected_thresholds)

    cluster_plot <- Seurat::DimPlot(
      seurat,
      cols = color_palette_cluster_seurat,
      label = TRUE,
      label.size = 3.5,
      reduction = reduction,
      group.by = "cluster_num"
    ) +
      theme(legend.position = "none") +
      xlab(
        paste0(
          ifelse(reduction == "tsne", "tSNE_1", "UMAP_1")
        )
      ) +
      ylab(
        paste0(
          ifelse(reduction == "tsne", "tSNE_2", "UMAP_2")
        )
      ) +
      labs(title = "Cluster") +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none"
      )

    combined_plots <- cluster_plot +
      auc_plots[[1]] +
      auc_plots[[2]] +
      auc_plots[[3]] +
      auc_plots[[4]] +
      auc_plots[[5]] +
      auc_plots[[6]] +
      auc_plots[[7]] +
      auc_plots[[8]] +
      plot_layout(
        ncol = 3,
        nrow = 3
      ) +
      plot_annotation(
        tag_levels = "a"
      )

    ggsave(
      file.path(
        output_dir,
        paste0("AUCell_combined_", reduction, ".pdf")
      ),
      combined_plots,
      width = 9.5,
      height = 10
    )

    ggsave(
      file.path(output_dir, paste0("AUCell_combined_", reduction, ".png")),
      combined_plots & theme(text = element_text(family = "Times New Roman")),
      width = 9.5,
      height = 10,
      dpi = 600
    )
    combined_plots_2 <- cluster_plot +
      auc_plots[[1]] +
      auc_plots[[2]] +
      auc_plots[[3]] +
      auc_plots[[4]] +
      auc_plots[[5]] +
      auc_plots[[6]] +
      auc_plots[[7]] +
      auc_plots[[8]] +
      plot_layout(
        ncol = 3,
        nrow = 3
      ) +
      plot_annotation(
        tag_levels = "a",
        tag_prefix = "(",
        tag_suffix = ")"
      )
    ggsave(
      file.path(output_dir, paste0("AUCell_combined_", reduction, "_2.pdf")),
      combined_plots_2,
      width = 9.5,
      height = 10
    )

    ggsave(
      file.path(output_dir, paste0("AUCell_combined_", reduction, "_2.png")),
      combined_plots_2 & theme(text = element_text(family = "Times New Roman")),
      width = 9.5,
      height = 10,
      dpi = 600
    )
  }

  auc_stats_data <- prepare_auc_data(
    seurat,
    selected_thresholds,
    color_palette_method
  )
  threshold_data <- threshold_proportion_data(
    seurat,
    selected_thresholds,
    color_palette_method
  )

  method_plots <- lapply(
    levels(auc_stats_data$GeneSet),
    function(gene_set_name) {
      method_name <- strsplit(gene_set_name, " ")[[1]][1]
      if (grepl("Random", method_name)) {
        method_name <- gene_set_name
      }

      subset_data <- auc_stats_data[auc_stats_data$GeneSet == gene_set_name, ]
      subset_data <- subset_data %>%
        group_by(CellType) %>%
        sample_frac(0.3) %>%
        ungroup()

      p_box <- ggplot(subset_data, aes(x = CellType, y = AUC)) +
        geom_boxplot(
          aes(fill = CellType),
          width = 0.8,
          outlier.shape = NA
        ) +
        geom_jitter(
          width = 0.25,
          size = 0.3,
          alpha = 0.4
        ) +
        scale_fill_manual(
          values = color_palette_cluster
        ) +
        theme_bw() +
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none",
          plot.margin = margin(-5, 5, -5, 5)
        ) +
        labs(
          x = NULL,
          y = "AUC"
        )

      method_data <- threshold_data[threshold_data$GeneSet == gene_set_name, ]
      method_data$Status <- factor(
        method_data$Status,
        levels = c("Above threshold", "Below threshold")
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

        above_data <- method_data[method_data$Status == "Above threshold", ]
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
          width = 0.8
        ) +
        scale_fill_manual(
          values = c(
            "Above threshold" = "black",
            "Below threshold" = "gray50"
          )
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
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          plot.margin = margin(-5, 5, 10, 5)
        ) +
        labs(x = "Cell type", y = "Proportion")

      combined_plot <- p_box / p_threshold +
        plot_layout(heights = c(2, 1)) &
        theme(panel.spacing = unit(0, "mm"))

      return(combined_plot)
    }
  )

  for (reduction in reductions) {
    auc_plots <- lapply(
      names(selected_thresholds), function(gene_set_name) {
        plot_auc_dimred(
          seurat,
          gene_set_name,
          reduction = reduction,
          show_labels = TRUE,
          point_size = 0.5,
          color_palette_method
        )
      }
    )
    names(auc_plots) <- names(selected_thresholds)

    final_plot_all <- auc_plots[[1]] +
      auc_plots[[2]] +
      auc_plots[[3]] +
      auc_plots[[4]] +
      method_plots[[1]] +
      method_plots[[2]] +
      method_plots[[3]] +
      method_plots[[4]] +
      auc_plots[[5]] +
      auc_plots[[6]] +
      auc_plots[[7]] +
      auc_plots[[8]] +
      method_plots[[5]] +
      method_plots[[6]] +
      method_plots[[7]] +
      method_plots[[8]] +
      plot_layout(
        guides = "collect",
        ncol = 4,
        nrow = 4
      ) & theme(
      legend.position = "bottom"
    )

    ggsave(
      file.path(output_dir, paste0("AUCell_stats_plots_", reduction, ".pdf")),
      final_plot_all,
      width = 12,
      height = 12
    )

    ggsave(
      file.path(output_dir, paste0("AUCell_stats_plots_", reduction, ".png")),
      final_plot_all & theme(text = element_text(family = "Times New Roman")),
      width = 12,
      height = 12,
      dpi = 600,
      bg = "white"
    )

    final_plot_all_2 <- auc_plots[[1]] +
      auc_plots[[2]] +
      auc_plots[[3]] +
      auc_plots[[4]] +
      method_plots[[1]] +
      method_plots[[2]] +
      method_plots[[3]] +
      method_plots[[4]] +
      auc_plots[[5]] +
      auc_plots[[6]] +
      auc_plots[[7]] +
      auc_plots[[8]] +
      method_plots[[5]] +
      method_plots[[6]] +
      method_plots[[7]] +
      method_plots[[8]] +
      plot_layout(
        guides = "collect",
        ncol = 4,
        nrow = 4
      )

    ggsave(
      file.path(output_dir, paste0("AUCell_stats_plots_2_", reduction, ".pdf")),
      final_plot_all_2,
      width = 12,
      height = 12
    )

    ggsave(
      file.path(output_dir, paste0("AUCell_stats_plots_2_", reduction, ".png")),
      final_plot_all_2 & theme(text = element_text(family = "Times New Roman")),
      width = 12,
      height = 12,
      dpi = 600,
      bg = "white"
    )
  }

  return(seurat)
}
