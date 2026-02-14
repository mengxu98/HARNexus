source("code/functions/prepare_env.R")

sample_pairs <- list(
  list(human = "h4", chimp = "c4"),
  list(human = "h3", chimp = "c2"),
  list(human = "h1", chimp = "c1")
)

human_color <- "#3271AE"
chimp_color <- "#D11A2D"

heatmap_color <- "#519673"

height_col <- 0.9
width_row <- 0.95

for (pair in sample_pairs) {
  human_sample <- pair$human
  chimp_sample <- pair$chimp

  log_message(
    "Similarity visualization for {.val {c(human_sample, chimp_sample)}}"
  )

  res_dir <- paste0(
    "results/species_networks/", human_sample, "_", chimp_sample, "/"
  )
  fig_dir <- check_dir(
    paste0("figures/species_networks/", human_sample, "_", chimp_sample)
  )

  summary_df <- read.csv(
    file.path(res_dir, "jaccard_similarity.csv"),
    stringsAsFactors = FALSE
  )

  similarity_human <- as.matrix(
    read.csv(
      file.path(res_dir, "similarity_matrix_human.csv"),
      row.names = 1,
      check.names = FALSE
    )
  )
  similarity_chimp <- as.matrix(
    read.csv(
      file.path(res_dir, "similarity_matrix_chimp.csv"),
      row.names = 1,
      check.names = FALSE
    )
  )
  similarity_vector_human_chimp <- read.csv(
    file.path(res_dir, "similarity_vector_human_chimp.csv"),
    stringsAsFactors = FALSE
  )

  celltypes <- rownames(similarity_human)
  n_ct <- length(celltypes)

  summary_df <- summary_df[match(celltypes, summary_df$CellType), ]
  similarity_vector_human_chimp <- similarity_vector_human_chimp[match(
    celltypes, similarity_vector_human_chimp$CellType
  ), ]
  celltypes_all <- similarity_vector_human_chimp$CellType
  similarity_vector_human_chimp <- similarity_vector_human_chimp$Jaccard_Edges_HumanChimp
  names(similarity_vector_human_chimp) <- celltypes_all
  chimp_celltypes_exist <- !is.na(summary_df$Chimp_Edges)

  combined_matrix <- matrix(
    NA,
    nrow = n_ct,
    ncol = n_ct,
    dimnames = list(celltypes, celltypes)
  )

  for (i in seq_len(n_ct)) {
    for (j in seq_len(n_ct)) {
      if (i < j) {
        combined_matrix[i, j] <- similarity_human[i, j]
      }
    }
  }

  for (i in seq_len(n_ct)) {
    for (j in seq_len(n_ct)) {
      if (i > j) {
        if (chimp_celltypes_exist[i] && chimp_celltypes_exist[j]) {
          combined_matrix[i, j] <- similarity_chimp[i, j]
        } else {
          combined_matrix[i, j] <- NA
        }
      }
    }
  }

  for (i in seq_len(n_ct)) {
    if (chimp_celltypes_exist[i]) {
      combined_matrix[i, i] <- similarity_vector_human_chimp[i]
    } else {
      combined_matrix[i, i] <- NA
    }
  }

  original_rownames <- rownames(combined_matrix)
  original_colnames <- colnames(combined_matrix)
  reversed_rownames <- rev(original_rownames)

  combined_matrix <- combined_matrix[rev(seq_len(nrow(combined_matrix))), , drop = FALSE]

  rownames(combined_matrix) <- reversed_rownames

  chimp_edges_vec <- ifelse(
    is.na(summary_df$Chimp_Edges), 0, summary_df$Chimp_Edges
  )
  chimp_tfs_vec <- ifelse(
    is.na(summary_df$Chimp_TFs), 0, summary_df$Chimp_TFs
  )
  chimp_genes_vec <- ifelse(
    is.na(summary_df$Chimp_Genes), 0, summary_df$Chimp_Genes
  )
  chimp_edges_vec <- chimp_edges_vec[match(original_colnames, summary_df$CellType)]
  chimp_tfs_vec <- chimp_tfs_vec[match(original_colnames, summary_df$CellType)]
  chimp_genes_vec <- chimp_genes_vec[match(original_colnames, summary_df$CellType)]

  human_edges_vec <- summary_df$Human_Edges[match(reversed_rownames, summary_df$CellType)]
  human_tfs_vec <- summary_df$Human_TFs[match(reversed_rownames, summary_df$CellType)]
  human_genes_vec <- summary_df$Human_Genes[match(reversed_rownames, summary_df$CellType)]

  jaccard_tfs_vec <- summary_df$Jaccard_TFs[match(original_colnames, summary_df$CellType)]
  jaccard_targets_vec <- summary_df$Jaccard_Targets[match(original_colnames, summary_df$CellType)]

  log_message("Creating Jaccard similarity heatmap...")

  chimp_edges_max <- max(chimp_edges_vec, na.rm = TRUE)
  chimp_tfs_max <- max(chimp_tfs_vec, na.rm = TRUE)
  chimp_genes_max <- max(chimp_genes_vec, na.rm = TRUE)
  chimp_edges_breaks <- pretty(c(0, chimp_edges_max), n = 3)
  chimp_tfs_breaks <- pretty(c(0, chimp_tfs_max), n = 3)
  chimp_genes_breaks <- pretty(c(0, chimp_genes_max), n = 3)
  chimp_edges_labels <- format(chimp_edges_breaks[-length(chimp_edges_breaks)], scientific = FALSE, big.mark = ",")
  chimp_tfs_labels <- format(chimp_tfs_breaks[-length(chimp_tfs_breaks)], scientific = FALSE, big.mark = ",")
  chimp_genes_labels <- format(chimp_genes_breaks[-length(chimp_genes_breaks)], scientific = FALSE, big.mark = ",")


  bar_colors_main <- sapply(original_colnames, function(ct) {
    if (ct %in% names(color_celltypes)) {
      color_celltypes[ct]
    } else {
      "#CCCCCC"
    }
  })

  celltype_values_top <- seq_along(original_colnames)
  col_mapping_top <- bar_colors_main
  names(col_mapping_top) <- as.character(celltype_values_top)

  top_color_bar <- ComplexHeatmap::columnAnnotation(
    "CellType" = ComplexHeatmap::anno_simple(
      celltype_values_top,
      col = col_mapping_top,
      height = unit(0.3, "cm")
    ),
    height = unit(0.3, "cm"),
    show_annotation_name = FALSE,
    annotation_name_gp = gpar(fontsize = 0)
  )

  col_anno <- ComplexHeatmap::columnAnnotation(
    "Chimpanzee edges" = ComplexHeatmap::anno_barplot(
      chimp_edges_vec,
      baseline = 0,
      bar_width = 0.5,
      gp = gpar(fill = chimp_color),
      height = unit(height_col, "cm"),
      axis_param = list(
        labels = chimp_edges_labels
      )
    ),
    "Chimpanzee TFs" = ComplexHeatmap::anno_barplot(
      chimp_tfs_vec,
      baseline = 0,
      bar_width = 0.5,
      gp = gpar(fill = chimp_color),
      height = unit(height_col, "cm"),
      axis_param = list(
        labels = chimp_tfs_labels
      )
    ),
    "Chimpanzee target genes" = ComplexHeatmap::anno_barplot(
      chimp_genes_vec,
      baseline = 0,
      bar_width = 0.5,
      gp = gpar(fill = chimp_color),
      height = unit(height_col, "cm"),
      axis_param = list(
        labels = chimp_genes_labels
      )
    ),
    height = unit(3 * height_col, "cm"),
    gap = unit(0.1, "cm"),
    annotation_name_gp = gpar(fontsize = 8)
  )

  bottom_color_bar <- ComplexHeatmap::columnAnnotation(
    "CellType" = ComplexHeatmap::anno_simple(
      celltype_values_top,
      col = col_mapping_top,
      height = unit(0.15, "cm")
    ),
    height = unit(0.15, "cm"),
    show_annotation_name = FALSE,
    annotation_name_gp = gpar(fontsize = 0)
  )

  human_edges_max <- max(human_edges_vec, na.rm = TRUE)
  human_tfs_max <- max(human_tfs_vec, na.rm = TRUE)
  human_genes_max <- max(human_genes_vec, na.rm = TRUE)
  human_edges_breaks <- pretty(c(0, human_edges_max), n = 3)
  human_tfs_breaks <- pretty(c(0, human_tfs_max), n = 3)
  human_genes_breaks <- pretty(c(0, human_genes_max), n = 3)
  human_edges_labels <- format(
    human_edges_breaks[-length(human_edges_breaks)],
    scientific = FALSE, big.mark = ","
  )
  human_tfs_labels <- format(
    human_tfs_breaks[-length(human_tfs_breaks)],
    scientific = FALSE, big.mark = ","
  )
  human_genes_labels <- format(
    human_genes_breaks[-length(human_genes_breaks)],
    scientific = FALSE, big.mark = ","
  )

  bar_colors_row <- sapply(reversed_rownames, function(ct) {
    if (ct %in% names(color_celltypes)) {
      color_celltypes[ct]
    } else {
      "#CCCCCC"
    }
  })

  celltype_values_row <- seq_along(reversed_rownames)
  col_mapping_row <- bar_colors_row
  names(col_mapping_row) <- as.character(celltype_values_row)

  row_anno <- ComplexHeatmap::rowAnnotation(
    "Human\ntarget\ngenes" = ComplexHeatmap::anno_barplot(
      human_genes_vec,
      baseline = 0,
      bar_width = 0.5,
      gp = gpar(fill = human_color),
      width = unit(width_row, "cm"),
      axis_param = list(
        labels = human_genes_labels
      )
    ),
    "Human\nTFs" = ComplexHeatmap::anno_barplot(
      human_tfs_vec,
      baseline = 0,
      bar_width = 0.5,
      gp = gpar(fill = human_color),
      width = unit(width_row, "cm"),
      axis_param = list(
        labels = human_tfs_labels
      )
    ),
    "Human\nedges" = ComplexHeatmap::anno_barplot(
      human_edges_vec,
      baseline = 0,
      bar_width = 0.5,
      gp = gpar(fill = human_color),
      width = unit(width_row, "cm"),
      axis_param = list(
        labels = human_edges_labels
      )
    ),
    " " = ComplexHeatmap::anno_simple(
      celltype_values_row,
      col = col_mapping_row,
      width = unit(0.15, "cm")
    ),
    width = unit(3 * width_row + 0.15, "cm"),
    gap = unit(0.1, "cm"),
    annotation_name_gp = gpar(fontsize = 8)
  )

  col_fun <- circlize::colorRamp2(
    c(0, 0.25),
    colorRampPalette(c("white", heatmap_color))(2)
  )

  p_heatmap <- ComplexHeatmap::Heatmap(
    combined_matrix,
    name = "Jaccard\nsimilarity",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    col = col_fun,
    na_col = "white",
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 9),
    show_column_names = TRUE,
    column_names_rot = 30,
    column_names_gp = gpar(fontsize = 9),
    column_names_side = "bottom",
    row_names_side = "right",
    cell_fun = function(j, i, x, y, width, height, fill) {
      val <- combined_matrix[i, j]
      n_rows <- nrow(combined_matrix)
      original_i <- n_rows - i + 1
      original_j <- j
      if (original_i < original_j) {
        border_col <- human_color
      } else if (original_i > original_j) {
        border_col <- chimp_color
      } else {
        border_col <- "gray80"
      }
      if (is.na(val)) {
        grid.rect(
          x = x, y = y, width = width, height = height,
          gp = gpar(fill = "white", col = border_col, lwd = 1)
        )
        grid.text(
          "NA", x, y,
          gp = gpar(fontsize = 8, col = "gray60", fontface = "italic")
        )
      } else {
        grid.rect(
          x = x, y = y, width = width, height = height,
          gp = gpar(fill = NA, col = border_col, lwd = 1)
        )
        text_col <- if (val > 0.1) "white" else "gray40"
        grid.text(
          sprintf("%.2f", val), x, y,
          gp = gpar(fontsize = 8, col = text_col)
        )
      }
    },
    width = unit(3.5, "cm"),
    height = unit(3.5, "cm"),
    heatmap_legend_param = list(
      title = "Jaccard\nsimilarity",
      title_gp = gpar(fontsize = 9),
      legend_height = unit(1, "cm"),
      at = c(0, 0.25),
      labels = c("0.00", "0.25")
    ),
    right_annotation = row_anno,
    top_annotation = col_anno,
    bottom_annotation = bottom_color_bar
  )

  pdf(
    file.path(fig_dir, "similarity_matrix.pdf"),
    width = 6,
    height = 5
  )
  ComplexHeatmap::draw(p_heatmap)
  dev.off()

  jaccard_max <- max(c(jaccard_tfs_vec, jaccard_targets_vec), na.rm = TRUE)
  if (is.na(jaccard_max) || jaccard_max == 0) {
    jaccard_max <- 1
  }
  jaccard_breaks <- pretty(c(0, jaccard_max), n = 3)
  jaccard_labels <- format(jaccard_breaks[-length(jaccard_breaks)], digits = 2)

  mean_tfs <- mean(jaccard_tfs_vec, na.rm = TRUE)
  mean_targets <- mean(jaccard_targets_vec, na.rm = TRUE)
  has_mean_tfs <- !is.nan(mean_tfs) && length(jaccard_tfs_vec) > 0
  has_mean_targets <- !is.nan(mean_targets) && length(jaccard_targets_vec) > 0

  axis_at_tfs <- sort(unique(c(jaccard_breaks[-length(jaccard_breaks)], mean_tfs)))
  axis_labels_tfs <- format(round(axis_at_tfs, 2), nsmall = 2)
  axis_labels_tfs[axis_at_tfs <= 1e-9] <- ""
  if (has_mean_tfs) {
    idx_mean <- which.min(abs(axis_at_tfs - mean_tfs))
    axis_labels_tfs[idx_mean] <- paste0("mean = ", axis_labels_tfs[idx_mean])
  }
  axis_at_targets <- sort(unique(c(jaccard_breaks[-length(jaccard_breaks)], mean_targets)))
  axis_labels_targets <- format(round(axis_at_targets, 2), nsmall = 2)
  axis_labels_targets[axis_at_targets <= 1e-9] <- ""
  if (has_mean_targets) {
    idx_mean <- which.min(abs(axis_at_targets - mean_targets))
    axis_labels_targets[idx_mean] <- paste0("mean = ", axis_labels_targets[idx_mean])
  }

  jaccard_tfs_plot <- ifelse(
    is.na(jaccard_tfs_vec), 0, jaccard_tfs_vec
  )
  jaccard_targets_plot <- ifelse(
    is.na(jaccard_targets_vec), 0, jaccard_targets_vec
  )

  jaccard_tfs_na_indices <- which(is.na(jaccard_tfs_vec))
  jaccard_targets_na_indices <- which(is.na(jaccard_targets_vec))

  bar_colors <- sapply(
    original_colnames, function(ct) {
      if (ct %in% names(color_celltypes)) {
        color_celltypes[ct]
      } else {
        "#CCCCCC"
      }
    }
  )

  dummy_matrix <- matrix(0, nrow = 2, ncol = length(original_colnames))
  rownames(dummy_matrix) <- c("TFs", "Genes")
  colnames(dummy_matrix) <- original_colnames

  jaccard_anno <- ComplexHeatmap::columnAnnotation(
    "Jaccard of TFs\n(Human vs Chimpanzee)" = ComplexHeatmap::anno_barplot(
      jaccard_tfs_plot,
      baseline = 0,
      bar_width = 0.5,
      gp = gpar(fill = bar_colors),
      height = unit(2, "cm"),
      axis_param = list(
        labels = axis_labels_tfs,
        at = axis_at_tfs
      ),
      ylim = c(0, jaccard_max)
    ),
    "Jaccard of target genes\n(Human vs Chimpanzee)" = ComplexHeatmap::anno_barplot(
      jaccard_targets_plot,
      baseline = 0,
      bar_width = 0.5,
      gp = gpar(fill = bar_colors),
      height = unit(2, "cm"),
      axis_param = list(
        labels = axis_labels_targets,
        at = axis_at_targets
      ),
      ylim = c(0, jaccard_max)
    ),
    height = unit(2, "cm"),
    gap = unit(0.2, "cm"),
    annotation_name_gp = gpar(fontsize = 9)
  )

  celltype_values <- seq_along(original_colnames)
  col_mapping <- bar_colors
  names(col_mapping) <- as.character(celltype_values)

  p_jaccard <- ComplexHeatmap::Heatmap(
    dummy_matrix,
    name = "Jaccard\nsimilarity",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_heatmap_legend = FALSE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_names_rot = 30,
    column_names_gp = gpar(fontsize = 9),
    width = unit(3, "cm"),
    height = unit(0.1, "cm"),
    top_annotation = jaccard_anno
  )

  pdf(
    file.path(fig_dir, "similarity_barplots.pdf"),
    width = 4,
    height = 3
  )
  ht_draw <- ComplexHeatmap::draw(p_jaccard)

  if (length(jaccard_tfs_na_indices) > 0) {
    ComplexHeatmap::decorate_annotation(
      "Jaccard of TFs\n(Human vs Chimpanzee)",
      {
        for (idx in jaccard_tfs_na_indices) {
          n_bars <- length(jaccard_tfs_vec)
          x_pos <- (idx - 0.5) / n_bars
          grid.text(
            "NA",
            x = unit(x_pos, "npc"),
            y = unit(0.5, "npc"),
            gp = gpar(fontsize = 8, col = "gray60", fontface = "italic")
          )
        }
      }
    )
  }

  if (has_mean_tfs) {
    ComplexHeatmap::decorate_annotation(
      "Jaccard of TFs\n(Human vs Chimpanzee)",
      {
        grid.lines(
          x = unit(c(0, 1), "npc"),
          y = unit(rep(mean_tfs, 2), "native"),
          gp = gpar(lty = 2, col = "black", lwd = 1)
        )
      }
    )
  }

  if (length(jaccard_targets_na_indices) > 0) {
    ComplexHeatmap::decorate_annotation(
      "Jaccard of target genes\n(Human vs Chimpanzee)",
      {
        for (idx in jaccard_targets_na_indices) {
          n_bars <- length(jaccard_targets_vec)
          x_pos <- (idx - 0.5) / n_bars
          grid.text(
            "NA",
            x = unit(x_pos, "npc"),
            y = unit(0.5, "npc"),
            gp = gpar(fontsize = 8, col = "gray60", fontface = "italic")
          )
        }
      }
    )
  }

  if (has_mean_targets) {
    ComplexHeatmap::decorate_annotation(
      "Jaccard of target genes\n(Human vs Chimpanzee)",
      {
        grid.lines(
          x = unit(c(0, 1), "npc"),
          y = unit(rep(mean_targets, 2), "native"),
          gp = gpar(lty = 2, col = "black", lwd = 1)
        )
      }
    )
  }

  dev.off()

  log_message(
    "Completed similarity plot for {.val {pair}}",
    message_type = "success"
  )
}
