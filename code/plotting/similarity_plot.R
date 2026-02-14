source("code/functions/prepare_env.R")

fig_dir <- check_dir("figures/networks/")

regions <- c("Prefrontal cortex", "Cerebral cortex")
exclude_tfs_in_genes <- TRUE

for (region in regions) {
  network_file <- file.path(
    "results/networks/analysis", paste0("network_", region, ".csv")
  )

  if (!file.exists(network_file)) {
    log_message(
      "Network file not found: {.val {network_file}}, skipping",
      message_type = "warning"
    )
    next
  }

  log_message("Processing region: {.val {region}}")

  network_data <- read.csv(network_file, stringsAsFactors = FALSE)
  if (exclude_tfs_in_genes) {
    network_data <- network_data[!network_data$Target %in% network_data$TF, ]
  }

  similarity_matrix_stage <- read.csv(
    file.path(
      "results/networks/analysis",
      paste0("similarity_matrix_Stage_", region, ".csv")
    ),
    row.names = 1,
    check.names = FALSE
  )
  similarity_matrix_celltype <- read.csv(
    file.path(
      "results/networks/analysis",
      paste0("similarity_matrix_CellType_", region, ".csv")
    ),
    row.names = 1,
    check.names = FALSE
  )

  similarity_matrices <- list(
    Stage = as.matrix(similarity_matrix_stage),
    CellType = as.matrix(similarity_matrix_celltype)
  )

  for (type in names(similarity_matrices)) {
    similarity_matrix <- similarity_matrices[[type]]
    if (type == "Stage") {
      subgroups <- sort(unique(network_data$CellType))
      subgroup_colors <- color_celltypes[subgroups]
      missing_colors <- is.na(subgroup_colors)
      if (any(missing_colors)) {
        subgroup_colors[missing_colors] <- "#CCCCCC"
      }
      names(subgroup_colors) <- subgroups

      sim_cols <- colnames(similarity_matrix)
      sim_rows <- rownames(similarity_matrix)

      edges_by_subgroup <- matrix(0,
        nrow = length(subgroups),
        ncol = length(sim_cols),
        dimnames = list(subgroups, sim_cols)
      )
      for (g in sim_cols) {
        for (subg in subgroups) {
          edges_by_subgroup[subg, g] <- sum(
            network_data[[type]] == g & network_data$CellType == subg
          )
        }
      }

      tfs_by_subgroup <- matrix(0,
        nrow = length(subgroups) + 1,
        ncol = length(sim_rows),
        dimnames = list(c("Shared", subgroups), sim_rows)
      )
      for (g in sim_rows) {
        group_data_all <- network_data[network_data[[type]] == g, ]
        all_tfs <- unique(group_data_all$TF)

        tfs_by_subg_list <- lapply(subgroups, function(subg) {
          subg_data <- group_data_all[group_data_all$CellType == subg, ]
          unique(subg_data$TF)
        })
        names(tfs_by_subg_list) <- subgroups

        tf_counts <- table(unlist(tfs_by_subg_list))
        shared_tfs <- names(tf_counts[tf_counts >= 2])
        tfs_by_subgroup["Shared", g] <- length(shared_tfs)

        for (subg in subgroups) {
          subg_tfs <- tfs_by_subg_list[[subg]]
          specific_tfs <- setdiff(subg_tfs, shared_tfs)
          tfs_by_subgroup[subg, g] <- length(specific_tfs)
        }
      }

      genes_by_subgroup <- matrix(0,
        nrow = length(subgroups) + 1,
        ncol = length(sim_rows),
        dimnames = list(c("Shared", subgroups), sim_rows)
      )
      for (g in sim_rows) {
        group_data_all <- network_data[network_data[[type]] == g, ]
        all_genes <- unique(group_data_all$Target)

        genes_by_subg_list <- lapply(subgroups, function(subg) {
          subg_data <- group_data_all[group_data_all$CellType == subg, ]
          unique(subg_data$Target)
        })
        names(genes_by_subg_list) <- subgroups

        gene_counts <- table(unlist(genes_by_subg_list))
        shared_genes <- names(gene_counts[gene_counts >= 2])
        genes_by_subgroup["Shared", g] <- length(shared_genes)

        for (subg in subgroups) {
          subg_genes <- genes_by_subg_list[[subg]]
          specific_genes <- setdiff(subg_genes, shared_genes)
          genes_by_subgroup[subg, g] <- length(specific_genes)
        }
      }
    } else if (type == "CellType") {
      subgroups <- unique(network_data$Stage)
      subgroup_numbers <- as.numeric(gsub("S", "", subgroups))
      subgroups <- subgroups[order(subgroup_numbers)]

      sim_cols <- colnames(similarity_matrix)
      sim_rows <- rownames(similarity_matrix)

      edges_by_subgroup <- matrix(0,
        nrow = length(subgroups),
        ncol = length(sim_cols),
        dimnames = list(subgroups, sim_cols)
      )
      for (g in sim_cols) {
        for (subg in subgroups) {
          edges_by_subgroup[subg, g] <- sum(
            network_data[[type]] == g & network_data$Stage == subg
          )
        }
      }

      tfs_by_subgroup <- matrix(0,
        nrow = length(subgroups) + 1,
        ncol = length(sim_rows),
        dimnames = list(c("Shared", subgroups), sim_rows)
      )
      for (g in sim_rows) {
        group_data_all <- network_data[network_data[[type]] == g, ]
        all_tfs <- unique(group_data_all$TF)
        tfs_by_subg_list <- lapply(subgroups, function(subg) {
          subg_data <- group_data_all[group_data_all$Stage == subg, ]
          unique(subg_data$TF)
        })
        names(tfs_by_subg_list) <- subgroups
        tf_counts <- table(unlist(tfs_by_subg_list))
        shared_tfs <- names(tf_counts[tf_counts >= 2])
        tfs_by_subgroup["Shared", g] <- length(shared_tfs)

        for (subg in subgroups) {
          subg_tfs <- tfs_by_subg_list[[subg]]
          specific_tfs <- setdiff(subg_tfs, shared_tfs)
          tfs_by_subgroup[subg, g] <- length(specific_tfs)
        }
      }
      genes_by_subgroup <- matrix(0,
        nrow = length(subgroups) + 1,
        ncol = length(sim_rows),
        dimnames = list(c("Shared", subgroups), sim_rows)
      )
      for (g in sim_rows) {
        group_data_all <- network_data[network_data[[type]] == g, ]
        all_genes <- unique(group_data_all$Target)

        genes_by_subg_list <- lapply(subgroups, function(subg) {
          subg_data <- group_data_all[group_data_all$Stage == subg, ]
          unique(subg_data$Target)
        })
        names(genes_by_subg_list) <- subgroups

        gene_counts <- table(unlist(genes_by_subg_list))
        shared_genes <- names(gene_counts[gene_counts >= 2])
        genes_by_subgroup["Shared", g] <- length(shared_genes)

        for (subg in subgroups) {
          subg_genes <- genes_by_subg_list[[subg]]
          specific_genes <- setdiff(subg_genes, shared_genes)
          genes_by_subgroup[subg, g] <- length(specific_genes)
        }
      }

      subgroup_colors <- color_stages[subgroups]
    }

    sim_cols <- colnames(similarity_matrix)
    sim_rows <- rownames(similarity_matrix)

    edges_by_subgroup <- edges_by_subgroup[, sim_cols, drop = FALSE]
    tfs_by_subgroup <- tfs_by_subgroup[, sim_rows, drop = FALSE]
    genes_by_subgroup <- genes_by_subgroup[, sim_rows, drop = FALSE]

    similarity_matrix[lower.tri(similarity_matrix, diag = FALSE)] <- NA

    tfs_subgroup_t <- t(tfs_by_subgroup)
    genes_subgroup_t <- t(genes_by_subgroup)
    tfs_subgroup_t <- tfs_subgroup_t[sim_rows, , drop = FALSE]
    genes_subgroup_t <- genes_subgroup_t[sim_rows, , drop = FALSE]
    expected_cols <- c("Shared", names(subgroup_colors))
    tfs_subgroup_t <- tfs_subgroup_t[, expected_cols, drop = FALSE]
    genes_subgroup_t <- genes_subgroup_t[, expected_cols, drop = FALSE]

    tfs_color_vec <- sapply(
      colnames(tfs_subgroup_t), function(col_name) {
        if (col_name == "Shared") {
          color_val <- "#CCCCCC"
        } else {
          color_val <- subgroup_colors[col_name]
          if (is.na(color_val) || is.null(color_val)) {
            color_val <- "#CCCCCC"
          }
        }
        color_val
      }
    )

    tfs_gp <- gpar(fill = tfs_color_vec)

    genes_color_vec <- sapply(
      colnames(genes_subgroup_t), function(col_name) {
        if (col_name == "Shared") {
          color_val <- "#CCCCCC"
        } else {
          color_val <- subgroup_colors[col_name]
          if (is.na(color_val) || is.null(color_val)) {
            color_val <- "#CCCCCC"
          }
        }
        color_val
      }
    )

    genes_gp <- gpar(fill = genes_color_vec)

    tfs_max <- max(rowSums(tfs_subgroup_t), na.rm = TRUE)
    genes_max <- max(rowSums(genes_subgroup_t), na.rm = TRUE)

    tfs_breaks <- pretty(c(0, tfs_max), n = 3)
    tfs_breaks <- tfs_breaks[-length(tfs_breaks)]
    tfs_labels <- format(tfs_breaks, scientific = FALSE, big.mark = ",")

    genes_breaks <- pretty(c(0, genes_max), n = 3)
    genes_breaks <- genes_breaks[-length(genes_breaks)]
    genes_labels <- format(genes_breaks, scientific = FALSE, big.mark = ",")

    row_anno <- ComplexHeatmap::rowAnnotation(
      "TFs count" = ComplexHeatmap::anno_barplot(
        tfs_subgroup_t,
        baseline = 0,
        gp = tfs_gp,
        width = unit(2, "cm"),
        axis_param = list(
          at = tfs_breaks,
          labels = tfs_labels
        )
      ),
      "Target genes count" = ComplexHeatmap::anno_barplot(
        genes_subgroup_t,
        baseline = 0,
        gp = genes_gp,
        width = unit(2, "cm"),
        axis_param = list(
          at = genes_breaks,
          labels = genes_labels
        )
      ),
      " " = ComplexHeatmap::anno_simple(
        seq_along(sim_rows),
        col = structure(
          sapply(sim_rows, function(row_name) {
            if (type == "Stage") {
              if (row_name %in% names(color_stages)) {
                color_stages[row_name]
              } else {
                "#CCCCCC"
              }
            } else {
              if (row_name %in% names(color_celltypes)) {
                color_celltypes[row_name]
              } else {
                "#CCCCCC"
              }
            }
          }),
          names = as.character(seq_along(sim_rows))
        ),
        width = unit(0.15, "cm")
      ),
      width = unit(4.5, "cm"),
      annotation_name_gp = gpar(fontsize = 8)
    )

    edges_subgroup_t <- t(edges_by_subgroup)
    edges_subgroup_t <- edges_subgroup_t[sim_cols, , drop = FALSE]
    edges_subgroup_t <- edges_subgroup_t[, names(subgroup_colors), drop = FALSE]
    edges_color_vec <- sapply(colnames(edges_subgroup_t), function(subg) {
      color_val <- subgroup_colors[subg]
      if (is.na(color_val) || is.null(color_val)) {
        color_val <- "#CCCCCC"
      }
      color_val
    })

    edges_gp <- gpar(fill = edges_color_vec)
    edges_max <- max(rowSums(edges_subgroup_t), na.rm = TRUE)
    edges_breaks <- pretty(c(0, edges_max), n = 3)
    edges_breaks <- edges_breaks[-length(edges_breaks)]
    edges_labels <- format(edges_breaks, scientific = FALSE, big.mark = ",")

    col_anno <- ComplexHeatmap::columnAnnotation(
      "Edges count" = ComplexHeatmap::anno_barplot(
        edges_subgroup_t,
        baseline = 0,
        gp = edges_gp,
        height = unit(2, "cm"),
        axis_param = list(
          at = edges_breaks,
          labels = edges_labels
        )
      ),
      " " = ComplexHeatmap::anno_simple(
        seq_along(sim_cols),
        col = structure(
          sapply(sim_cols, function(col_name) {
            if (type == "Stage") {
              if (col_name %in% names(color_stages)) {
                color_stages[col_name]
              } else {
                "#CCCCCC"
              }
            } else {
              if (col_name %in% names(color_celltypes)) {
                color_celltypes[col_name]
              } else {
                "#CCCCCC"
              }
            }
          }),
          names = as.character(seq_along(sim_cols))
        ),
        height = unit(0.15, "cm")
      ),
      height = unit(2.5, "cm"),
      annotation_name_gp = gpar(fontsize = 8)
    )

    if (type == "Stage") {
      col_colors <- sapply(sim_cols, function(stage) {
        if (stage %in% names(color_stages)) {
          color_stages[stage]
        } else {
          "#CCCCCC"
        }
      })
    } else if (type == "CellType") {
      col_colors <- sapply(sim_cols, function(ct) {
        if (ct %in% names(color_celltypes)) {
          color_celltypes[ct]
        } else {
          "#CCCCCC"
        }
      })
    }

    celltype_values_col <- seq_along(sim_cols)
    col_mapping_col <- col_colors
    names(col_mapping_col) <- as.character(celltype_values_col)

    col_fun <- circlize::colorRamp2(
      c(0, 0.25, 0.5, 0.75, 1),
      colorRampPalette(c("white", "#2177B8"))(5)
    )

    if (type == "Stage") {
      size <- 7
      w <- 9
      h <- 5
    } else if (type == "CellType") {
      size <- 5
      w <- 10
      h <- 4
    }

    p_sim <- ComplexHeatmap::Heatmap(
      similarity_matrix,
      name = "Jaccard\nsimilarity",
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      col = col_fun,
      show_column_names = FALSE,
      na_col = "white",
      cell_fun = function(j, i, x, y, width, height, fill) {
        val <- similarity_matrix[i, j]
        if (is.na(val)) {
          grid.rect(
            x = x, y = y, width = width, height = height,
            gp = gpar(fill = "white", col = NA, lwd = 0)
          )
        } else {
          grid.rect(
            x = x, y = y, width = width, height = height,
            gp = gpar(fill = NA, col = "gray80", lwd = 0.5)
          )
          text_col <- if (val > 0.5) "white" else "gray40"
          grid.text(
            sprintf("%.2f", val), x, y,
            gp = gpar(fontsize = 8, col = text_col)
          )
        }
      },
      width = unit(size, "cm"),
      height = unit(size, "cm"),
      column_title = thisutils::capitalize(
        type,
        force_tolower = TRUE
      ),
      heatmap_legend_param = list(
        title = "Jaccard\nsimilarity",
        title_gp = gpar(fontsize = 10),
        legend_height = unit(4, "cm"),
        at = c(0, 0.25, 0.5, 0.75, 1),
        labels = c("0.00", "0.25", "0.50", "0.75", "1.00")
      ),
      right_annotation = row_anno,
      top_annotation = col_anno
    )

    legend_list <- list(
      ComplexHeatmap::Legend(
        title = "Edges/TFs/Target genes type",
        title_gp = gpar(fontsize = 10),
        at = names(subgroup_colors),
        labels = names(subgroup_colors),
        legend_gp = gpar(fill = subgroup_colors),
        ncol = if (length(subgroup_colors) > 9) 2 else 1
      ),
      ComplexHeatmap::Legend(
        title = "Shared TFs/Target genes type",
        title_gp = gpar(fontsize = 10),
        at = c("Shared"),
        labels = c("Shared"),
        legend_gp = gpar(fill = c("#CCCCCC")),
        ncol = 1
      )
    )

    pdf(
      file.path(
        fig_dir,
        paste0("similarity_matrix_", type, "_", region, ".pdf")
      ),
      width = w, height = h
    )
    ComplexHeatmap::draw(p_sim, annotation_legend_list = legend_list)
    dev.off()

    calculate_jaccard <- function(set1, set2) {
      if (length(set1) == 0 && length(set2) == 0) {
        return(0.0)
      }
      intersection <- length(intersect(set1, set2))
      union_len <- length(union(set1, set2))
      return(ifelse(union_len > 0, intersection / union_len, 0.0))
    }

    tfs_sets <- list()
    genes_sets <- list()
    for (g in sim_rows) {
      group_data_all <- network_data[network_data[[type]] == g, ]
      tfs_sets[[g]] <- unique(group_data_all$TF)
      genes_sets[[g]] <- unique(group_data_all$Target)
    }

    tfs_jaccard_vec <- numeric(length(sim_rows))
    genes_jaccard_vec <- numeric(length(sim_rows))
    names(tfs_jaccard_vec) <- sim_rows
    names(genes_jaccard_vec) <- sim_rows

    for (i in seq_along(sim_rows)) {
      g_i <- sim_rows[i]
      tfs_i <- tfs_sets[[g_i]]
      genes_i <- genes_sets[[g_i]]

      tfs_jaccards <- numeric(0)
      genes_jaccards <- numeric(0)
      for (j in seq_along(sim_rows)) {
        if (i != j) {
          g_j <- sim_rows[j]
          tfs_j <- tfs_sets[[g_j]]
          genes_j <- genes_sets[[g_j]]
          tfs_jaccards <- c(tfs_jaccards, calculate_jaccard(tfs_i, tfs_j))
          genes_jaccards <- c(
            genes_jaccards, calculate_jaccard(genes_i, genes_j)
          )
        }
      }
      tfs_jaccard_vec[i] <- mean(tfs_jaccards, na.rm = TRUE)
      genes_jaccard_vec[i] <- mean(genes_jaccards, na.rm = TRUE)
    }

    tfs_breaks <- pretty(c(0, 1), n = 3)
    genes_breaks <- pretty(c(0, 1), n = 3)
    tfs_labels <- format(tfs_breaks[-length(tfs_breaks)], digits = 2)
    genes_labels <- format(genes_breaks[-length(genes_breaks)], digits = 2)

    mean_tfs <- mean(tfs_jaccard_vec, na.rm = TRUE)
    mean_genes <- mean(genes_jaccard_vec, na.rm = TRUE)
    has_mean_tfs <- !is.nan(mean_tfs) && length(tfs_jaccard_vec) > 0
    has_mean_genes <- !is.nan(mean_genes) && length(genes_jaccard_vec) > 0

    axis_at_tfs <- sort(unique(c(tfs_breaks[-length(tfs_breaks)], mean_tfs)))
    axis_labels_tfs <- format(round(axis_at_tfs, 2), nsmall = 2)
    axis_labels_tfs[axis_at_tfs <= 1e-9] <- ""
    if (has_mean_tfs) {
      idx_mean <- which.min(abs(axis_at_tfs - mean_tfs))
      axis_labels_tfs[idx_mean] <- paste0("mean = ", axis_labels_tfs[idx_mean])
    }
    axis_at_genes <- sort(unique(c(genes_breaks[-length(genes_breaks)], mean_genes)))
    axis_labels_genes <- format(round(axis_at_genes, 2), nsmall = 2)
    axis_labels_genes[axis_at_genes <= 1e-9] <- ""
    if (has_mean_genes) {
      idx_mean <- which.min(abs(axis_at_genes - mean_genes))
      axis_labels_genes[idx_mean] <- paste0("mean = ", axis_labels_genes[idx_mean])
    }

    if (type == "Stage") {
      bar_colors <- sapply(sim_rows, function(stage) {
        if (stage %in% names(color_stages)) {
          color_stages[stage]
        } else {
          "#CCCCCC"
        }
      })
    } else if (type == "CellType") {
      bar_colors <- sapply(sim_rows, function(ct) {
        if (ct %in% names(color_celltypes)) {
          color_celltypes[ct]
        } else {
          "#CCCCCC"
        }
      })
    }

    dummy_matrix <- matrix(0, nrow = 2, ncol = length(sim_rows))
    rownames(dummy_matrix) <- c("TFs", "Genes")
    colnames(dummy_matrix) <- sim_rows

    barplots_anno <- ComplexHeatmap::columnAnnotation(
      "Jaccard of TFs" = ComplexHeatmap::anno_barplot(
        tfs_jaccard_vec,
        baseline = 0,
        bar_width = 0.5,
        gp = gpar(fill = bar_colors),
        height = unit(2, "cm"),
        axis_param = list(
          labels = axis_labels_tfs,
          at = axis_at_tfs
        ),
        ylim = c(0, 1)
      ),
      "Jaccard of target genes" = ComplexHeatmap::anno_barplot(
        genes_jaccard_vec,
        baseline = 0,
        bar_width = 0.5,
        gp = gpar(fill = bar_colors),
        height = unit(2, "cm"),
        axis_param = list(
          labels = axis_labels_genes,
          at = axis_at_genes
        ),
        ylim = c(0, 1)
      ),
      height = unit(2, "cm"),
      gap = unit(0.2, "cm"),
      annotation_name_gp = gpar(fontsize = 7)
    )

    p_barplots <- ComplexHeatmap::Heatmap(
      dummy_matrix,
      name = "Count",
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_heatmap_legend = FALSE,
      show_row_names = FALSE,
      show_column_names = TRUE,
      column_names_rot = 45,
      column_names_gp = gpar(fontsize = 8),
      width = unit(4, "cm"),
      height = unit(0.2, "cm"),
      top_annotation = barplots_anno
    )

    pdf(
      file.path(
        fig_dir,
        paste0("similarity_barplots_", type, "_", region, ".pdf")
      ),
      width = 4,
      height = 3
    )
    ComplexHeatmap::draw(p_barplots)

    if (has_mean_tfs) {
      ComplexHeatmap::decorate_annotation(
        "Jaccard of TFs",
        {
          grid.lines(
            x = unit(c(0, 1), "npc"),
            y = unit(rep(mean_tfs, 2), "native"),
            gp = gpar(lty = 2, col = "black", lwd = 1)
          )
        }
      )
    }
    if (has_mean_genes) {
      ComplexHeatmap::decorate_annotation(
        "Jaccard of target genes",
        {
          grid.lines(
            x = unit(c(0, 1), "npc"),
            y = unit(rep(mean_genes, 2), "native"),
            gp = gpar(lty = 2, col = "black", lwd = 1)
          )
        }
      )
    }

    dev.off()

    log_message("Saved similarity plot")
  }
}
