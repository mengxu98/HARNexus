source("code/functions/prepare_env.R")

res_dir <- check_dir("results/pfc_astrocytes")
fig_dir <- check_dir(
  "figures/pfc_astrocytes"
)

selected_stages <- c("S6", "S7", "S11", "S12")
stage_colors <- color_stages[selected_stages]


subcelltype_colors <- c(
  "Early astrocytes" = "#EEB8C3",
  "Middle astrocytes" = "#EE4866",
  "Late astrocytes" = "#B81A35"
)

feature_cluster_colors <- c(
  "C1" = "#0AA344",
  "C2" = "#F9BD10"
)

# Obtain marker genes from doi: https://doi.org/10.1038/s41556-024-01583-9
marker_genes_split <- data.frame(
  Celltype = c(
    rep("Early astrocytes", 3),
    rep("Middle astrocytes", 3),
    rep("Late astrocytes", 3)
  ),
  Genes = c(
    # Early / progenitor-like
    "MKI67", "TOP2A", "NES",
    # Middle / immature
    "POU3F3", "SOX2", "LHX2",
    # Late / mature
    "AQP4", "GJA1", "GRM3"
  )
)
objects_target <- readRDS(
  file.path(res_dir, "object_pfc_astrocytes_target.rds")
)

for (objects_name in c("var", "target")) {
  log_message("Processing {.val {objects_name}}...")

  title_name <- if (objects_name == "var") "Variable genes" else "Target genes"
  objects <- readRDS(
    file.path(res_dir, paste0("object_pfc_astrocytes_", objects_name, ".rds"))
  )
  objects$Stage <- factor(
    objects$Stage,
    levels = selected_stages
  )

  p_list <- list()
  for (by in c("seurat_clusters", "Sample", "Age")) {
    log_message("Plotting {.val {by}}...")
    p_list[[by]] <- CellDimPlot(
      objects,
      reduction = "umap",
      group.by = by,
      xlab = "UMAP_1",
      ylab = "UMAP_2"
    )
  }
  p1 <- p_list$Sample + p_list$Age + p_list$seurat_clusters
  ggsave(
    file.path(fig_dir, paste0("cell_dim_plot_", objects_name, ".pdf")),
    p1,
    width = 14,
    height = 4
  )

  p2 <- CellDimPlot(
    objects,
    group.by = "Stage",
    reduction = "umap",
    pt.size = 3,
    raster = TRUE,
    palcolor = stage_colors,
    title = title_name,
    xlab = "UMAP_1",
    ylab = "UMAP_2"
  )
  p3 <- CellDimPlot(
    objects,
    group.by = "Celltype",
    reduction = "umap",
    pt.size = 3,
    raster = TRUE,
    palcolor = subcelltype_colors,
    xlab = "UMAP_1",
    ylab = "UMAP_2"
  )
  p4 <- p2 + p3
  ggsave(
    file.path(
      fig_dir,
      paste0("cell_dim_plot_celltype_", objects_name, ".pdf")
    ),
    p4,
    width = 8.5,
    height = 3
  )

  p5 <- FeatureDimPlot(
    objects,
    features = marker_genes_split$Genes,
    reduction = "umap",
    pt.size = 7,
    raster = TRUE,
    xlab = "UMAP_1",
    ylab = "UMAP_2",
    ncol = 3
  )
  ggsave(
    file.path(fig_dir, paste0("feature_dim_plot_", objects_name, ".pdf")),
    p5,
    width = 7, height = 7
  )

  GroupHeatmap(
    objects,
    features = marker_genes_split$Genes,
    group.by = "seurat_clusters",
    heatmap_palette = "Spectral",
    height = 2.5,
    width = 2.2,
    add_dot = TRUE,
    dot_size = unit(6, "mm"),
    nlabel = 0,
    show_row_names = TRUE,
    show_column_names = TRUE
  )
  ggsave(
    file.path(fig_dir, paste0("group_heatmap_", objects_name, ".pdf")),
    width = 6, height = 4
  )

  GroupHeatmap(
    objects,
    features = marker_genes_split$Genes,
    feature_split = marker_genes_split$Celltype,
    group.by = "Celltype",
    group_palcolor = subcelltype_colors,
    feature_split_palcolor = subcelltype_colors,
    heatmap_palette = "Spectral",
    show_row_names = TRUE,
    nlabel = 0,
    height = 2.5,
    width = 1,
    add_dot = TRUE,
    dot_size = unit(8, "mm"),
  )
  ggsave(
    file.path(fig_dir, paste0("group_heatmap_celltype_", objects_name, ".pdf")),
    width = 7, height = 4
  )

  p6 <- FeatureDimPlot(
    objects,
    features = "Lineage1",
    title = title_name,
    palette = "viridis",
    # palcolor = c("#0AA344", "#F9BD10", "#D70440"),
    reduction = "umap",
    pt.size = 3,
    raster = TRUE,
    xlab = "UMAP_1",
    ylab = "UMAP_2"
  )
  ggsave(
    file.path(fig_dir, paste0("pseudotime_dim_plot_", objects_name, ".pdf")),
    p6,
    width = 4, height = 3.5
  )

  p7 <- CellDimPlot(
    objects,
    palcolor = stage_colors,
    group.by = "Stage",
    reduction = "umap",
    lineages = "Lineage1",
    pt.size = 3,
    raster = TRUE,
    title = title_name,
    lineages_span = 0.1,
    lineages_palcolor = "#F9BD10",
    xlab = "UMAP_1",
    ylab = "UMAP_2"
  )
  ggsave(
    file.path(fig_dir, paste0("slingshot_plot_", objects_name, ".pdf")),
    p7,
    width = 6,
    height = 3.5
  )
  p8 <- CellDimPlot(
    objects,
    palcolor = subcelltype_colors,
    group.by = "Celltype",
    reduction = "umap",
    lineages = "Lineage1",
    pt.size = 3,
    raster = TRUE,
    title = title_name,
    lineages_span = 0.1,
    lineages_palcolor = "#F9BD10",
    xlab = "UMAP_1",
    ylab = "UMAP_2"
  )
  ggsave(
    file.path(
      fig_dir,
      paste0("slingshot_plot_celltype_", objects_name, ".pdf")
    ),
    p8,
    width = 6,
    height = 3.5
  )

  file_name <- file.path(
    res_dir,
    paste0("dynamic_heatmap_", objects_name, ".rds")
  )
  if (!file.exists(file_name)) {
    features_file <- file.path(
      res_dir,
      paste0("dynamic_heatmap_", objects_name, "_feature.csv")
    )
    cluster_filter <- if (objects_name == "var") "C2" else "C1"
    features_from_csv <- if (file.exists(features_file)) {
      features_df <- read.csv(features_file, stringsAsFactors = FALSE)
      col2_name <- names(features_df)[2L]
      subset_df <- features_df[features_df[[col2_name]] == cluster_filter, , drop = FALSE]
      if (nrow(subset_df) > 0) {
        data.frame(
          gene = as.character(subset_df[[1L]]),
          group = cluster_filter,
          stringsAsFactors = FALSE
        )
      } else {
        NULL
      }
    } else {
      NULL
    }

    csv_dir <- "results/networks/har_csn_atlas/csv"
    region_name <- "Prefrontal cortex"
    celltype_name <- "Astrocytes"
    network_list <- lapply(
      selected_stages,
      function(s) {
        f <- file.path(
          csv_dir,
          paste0(region_name, "_", s, "_", celltype_name, ".csv")
        )
        if (!file.exists(f)) {
          return(NULL)
        }
        net <- read.csv(f, stringsAsFactors = FALSE)
        regulators <- unique(net$regulator)
        net <- net[!(net$target %in% regulators), ]
        net$Stage <- s
        net
      }
    )
    network_list <- network_list[!sapply(network_list, is.null)]
    network_data <- if (length(network_list) > 0) {
      do.call(rbind, network_list)
    } else {
      data.frame(
        regulator = character(),
        target = character(),
        weight = numeric()
      )
    }
    genes_csv <- features_from_csv[, 1L]
    network_data <- network_data[network_data$target %in% genes_csv, , drop = FALSE]
    network_data$abs_weight <- abs(network_data$weight)
    network_ordered <- network_data[order(-network_data$abs_weight), ]
    top10_edges <- head(network_ordered, 10)
    features_selected <- unique(top10_edges$target)
    write.csv(
      features_selected,
      file.path(
        res_dir,
        paste0("dynamic_heatmap_", objects_name, "_feature_selected.csv")
      ),
      quote = FALSE
    )

    ht <- DynamicHeatmap(
      objects,
      lineages = "Lineage1",
      features_label = features_selected,
      cell_annotation = c("Stage", "Celltype"),
      terms_fontsize = 10,
      heatmap_palette = "viridis",
      pseudotime_palette = "plasma",
      feature_split_palcolor = feature_cluster_colors,
      cell_annotation_palcolor = list(
        Stage = stage_colors,
        Celltype = subcelltype_colors
      ),
      anno_terms = TRUE,
      nlabel = 0,
      n_split = 2,
      width = 2.5,
      height = 4,
      use_raster = TRUE
    )
    saveRDS(
      ht,
      file_name
    )
  } else {
    ht <- readRDS(
      file_name
    )
  }
  ggsave(
    file.path(
      fig_dir,
      paste0("dynamic_heatmap_", objects_name, ".pdf")
    ),
    ht$plot,
    width = 10,
    height = 6
  )
  write.csv(
    ht$feature_split,
    file.path(
      res_dir,
      paste0("dynamic_heatmap_", objects_name, "_feature.csv")
    ),
    quote = FALSE
  )
}
