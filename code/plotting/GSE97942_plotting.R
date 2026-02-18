source("code/functions/prepare_env.R")

fig_dir <- check_dir("figures/gse97942")
object <- readRDS(
  "../../data/BrainData/processed/GSE97942/GSE97942_cerebellum_processed.rds"
)
celltype_levels <- sort(as.character(unique(object$CellType)))
object$Celltype <- factor(
  object$CellType,
  levels = celltype_levels
)

p1 <- CellDimPlot(
  object,
  reduction = "umap",
  group.by = "Sample",
  palcolor = colors9,
  pt.size = 3,
  raster = TRUE,
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)

p2 <- CellDimPlot(
  object,
  reduction = "umap",
  group.by = "seurat_clusters",
  pt.size = 3,
  raster = TRUE,
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)

p3 <- CellDimPlot(
  object,
  reduction = "umap",
  group.by = "Celltype",
  pt.size = 3,
  raster = TRUE,
  palcolor = color_celltypes[celltype_levels],
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)
p3 <- p1 / p2 / p3
ggsave(
  file.path(fig_dir, "dim_plots_celltype.pdf"),
  p3,
  width = 6,
  height = 10
)

marker_genes_split <- data.frame(
  Celltype = c(
    rep("Microglia", 3),
    rep("Oligodendrocytes", 3),
    rep("Oligodendrocyte progenitor cells", 5),
    rep("Inhibitory neurons", 3),
    rep("Astrocytes", 5),
    rep("Excitatory neurons", 1)
  ),
  Genes = c(
    "CX3CR1", "P2RY12", "CSF1R",
    "MOG", "MAG", "CLDN11",
    "PDGFRA", "CSPG4", "OLIG1", "OLIG2", "SOX10",
    "GAD1", "GAD2", "SLC6A1",
    "GFAP", "AQP4", "ALDH1L1", "FGFR3", "GJA1",
    "SLC17A7"
  )
)

marker_genes_split$Celltype <- factor(
  marker_genes_split$Celltype,
  levels = celltype_levels
)

celltype_colors <- color_celltypes[celltype_levels]
p4 <- FeatureDimPlot(
  object,
  features = marker_genes_split$Genes,
  reduction = "umap",
  ncol = 5,
  pt.size = 7,
  raster = TRUE,
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)
ggsave(
  file.path(fig_dir, "feature_dim_plot_markergenes.pdf"),
  p4,
  width = 13,
  height = 9
)

gh <- GroupHeatmap(
  object,
  exp_legend_title = "Z-score",
  features = marker_genes_split$Genes,
  feature_split = marker_genes_split$Celltype,
  group.by = "Celltype",
  group_palcolor = celltype_colors,
  cell_annotation_palcolor = celltype_colors,
  feature_split_palcolor = celltype_colors,
  heatmap_palette = "Spectral",
  height = 3.8,
  width = 2,
  add_dot = TRUE,
  dot_size = unit(12, "mm"),
  nlabel = 0,
  show_row_names = TRUE,
  ht_params = list(
    row_names_gp = gpar(fontface = "italic")
  )
)
pdf(
  file.path(fig_dir, "group_heatmap_markergenes.pdf"),
  width = 10.5,
  height = 4.6
)
print(gh$plot)
dev.off()
