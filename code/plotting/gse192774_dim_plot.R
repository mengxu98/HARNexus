source("code/functions/prepare_env.R")

sample_pairs <- list(
  list(human = "h4", chimp = "c4"),
  list(human = "h3", chimp = "c2"),
  list(human = "h1", chimp = "c1")
)

human_color <- "#3271AE"
chimp_color <- "#D11A2D"

species_colors <- c(
  "Human" = human_color, "Chimpanzee" = chimp_color
)

for (pair in sample_pairs) {
  human_sample <- pair$human
  chimp_sample <- pair$chimp

  log_message(
    "Loading data for {.val {pair}}..."
  )
  res_dir <- paste0(
    "results/species_networks/", human_sample, "_", chimp_sample, "/"
  )
  fig_dir <- check_dir(
    paste0(
      "figures/species_networks/", human_sample, "_", chimp_sample
    )
  )

  file_path <- file.path(res_dir, "merged_object.rds")
  if (!file.exists(file_path)) {
    human_obj <- readRDS(
      file.path(
        res_dir,
        paste0("human_", human_sample, "_object.rds")
      )
    )
    chimp_obj <- readRDS(
      file.path(
        res_dir,
        paste0("chimp_", chimp_sample, "_object.rds")
      )
    )
    human_obj$Species <- "Human"
    chimp_obj$Species <- "Chimpanzee"
    merged_obj <- merge(human_obj, y = chimp_obj)
    merged_obj <- NormalizeData(merged_obj)
    merged_obj <- FindVariableFeatures(merged_obj)
    merged_obj <- ScaleData(merged_obj)
    merged_obj <- RunPCA(merged_obj)
    merged_obj[["RNA"]] <- split(merged_obj[["RNA"]], merged_obj$Species)
    merged_obj <- IntegrateLayers(
      merged_obj,
      method = RPCAIntegration,
      orig.reduction = "pca",
      new.reduction = "integrated.rpca"
    )
    merged_obj <- RunUMAP(
      merged_obj,
      reduction = "integrated.rpca",
      dims = 1:30
    )
    saveRDS(merged_obj, file_path)
  } else {
    merged_obj <- readRDS(file_path)
  }

  p_celltype <- CellDimPlot(
    merged_obj,
    reduction = "umap",
    group.by = "CellType",
    palcolor = color_celltypes,
    label = FALSE,
    raster = TRUE,
    pt.size = 3,
    xlab = "UMAP_1",
    ylab = "UMAP_2"
  )

  p_species <- CellDimPlot(
    merged_obj,
    reduction = "umap",
    group.by = "Species",
    palcolor = species_colors,
    label = FALSE,
    raster = TRUE,
    pt.size = 3,
    xlab = "UMAP_1",
    ylab = "UMAP_2"
  )

  p_merged <- p_celltype + p_species

  ggsave(
    file.path(fig_dir, "dim_plots_merged.pdf"),
    p_merged,
    width = 9.5, height = 3
  )
}
