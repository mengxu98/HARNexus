rm(list = ls())
gc()

source("code/functions/prepare_env.R")

data_dir <- "../../data/BrainData/integration/"
fig_dir <- check_dir("results/figures/datasets")

objects <- readRDS(file.path(data_dir, "objects_celltype_plot.rds"))

p1 <- DimPlot(
  objects,
  reduction = "umap.unintegrated",
  group.by = "Dataset",
  cols = colors32,
  raster = TRUE
) +
  coord_fixed()
p2 <- DimPlot(
  objects,
  reduction = "umap.rpca",
  group.by = "Dataset",
  cols = colors32,
  raster = TRUE
) +
  coord_fixed()
p3 <- DimPlot(
  objects,
  reduction = "umap.harmony",
  group.by = "Dataset",
  cols = colors32,
  raster = TRUE
) +
  coord_fixed()
p4 <- p1 + p2 + p3 +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")
ggsave(
  file.path(fig_dir, "01_dim_plots_datasets.pdf"),
  p4,
  width = 22,
  height = 7
)

p5 <- DotPlot(
  objects,
  features = marker_genes,
  group.by = "seurat_clusters",
  cols = colors2,
  dot.scale = 5
) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.key.width = unit(0.2, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.text = element_text(size = 12)
  ) +
  coord_fixed()
ggsave(
  file.path(fig_dir, "02_dot_plot_markergenes.pdf"),
  p5,
  width = 9,
  height = 20
)

p6 <- FeaturePlot(
  objects,
  features = marker_genes,
  reduction = "umap.unintegrated",
  cols = colors2,
  ncol = 6,
  raster = TRUE
)
ggsave(
  file.path(fig_dir, "03_feature_plots_unintegrated.pdf"),
  p6,
  width = 20,
  height = 15
)
p7 <- FeaturePlot(
  objects,
  features = marker_genes,
  reduction = "umap.rpca",
  cols = colors2,
  ncol = 6,
  raster = TRUE
)
ggsave(
  file.path(fig_dir, "03_feature_plots_rpca.pdf"),
  p7,
  width = 20,
  height = 15
)
p8 <- FeaturePlot(
  objects,
  features = marker_genes,
  reduction = "umap.harmony",
  cols = colors2,
  ncol = 6,
  raster = TRUE
)
ggsave(
  file.path(fig_dir, "03_feature_plots_harmony.pdf"),
  p8,
  width = 20,
  height = 15
)

Idents(objects) <- "CellType"

Idents(objects) <- factor(
  Idents(objects),
  levels = names(color_celltype)
)

p8 <- DimPlot(
  objects,
  reduction = "umap.unintegrated",
  cols = colors128,
  group.by = "seurat_clusters",
  label = FALSE,
  raster = TRUE
) +
  coord_fixed()
p9 <- DimPlot(
  objects,
  reduction = "umap.unintegrated",
  cols = color_celltype,
  group.by = "CellType",
  label = FALSE,
  raster = TRUE
) +
  coord_fixed()
p10 <- DimPlot(
  objects,
  reduction = "umap.rpca",
  cols = colors128,
  group.by = "seurat_clusters",
  label = FALSE,
  raster = TRUE
) +
  coord_fixed()
p11 <- DimPlot(
  objects,
  reduction = "umap.rpca",
  cols = color_celltype,
  group.by = "CellType",
  label = FALSE,
  raster = TRUE
) +
  coord_fixed()
p12 <- DimPlot(
  objects,
  reduction = "umap.harmony",
  cols = colors128,
  group.by = "seurat_clusters",
  label = FALSE,
  raster = TRUE
) +
  coord_fixed()
p13 <- DimPlot(
  objects,
  reduction = "umap.harmony",
  cols = color_celltype,
  group.by = "CellType",
  label = FALSE,
  raster = TRUE
) +
  coord_fixed()

p14 <- p8 + p9 + plot_annotation(tag_levels = "A")
ggsave(
  file.path(fig_dir, "05_celltype_unintegrated.pdf"),
  p14,
  width = 18, height = 6
)
p15 <- p10 + p11 + plot_annotation(tag_levels = "A")
ggsave(
  file.path(fig_dir, "05_celltype_rpca.pdf"),
  p15,
  width = 18, height = 6
)
p16 <- p12 + p13 + plot_annotation(tag_levels = "A")
ggsave(
  file.path(fig_dir, "05_celltype_harmony.pdf"),
  p16,
  width = 18, height = 6
)


p9 <- DotPlot(
  objects,
  features = marker_genes,
  group.by = "CellType",
  dot.scale = 6,
  cols = colors2
) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.key.width = unit(0.2, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.text = element_text(size = 12)
  ) +
  labs(x = "Marker genes", y = "Cell type") +
  coord_fixed()

ggsave(
  file.path(fig_dir, "04_dot_plot_celltype.pdf"),
  p9,
  width = 11, height = 3
)
