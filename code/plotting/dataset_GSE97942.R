rm(list = ls())
gc()

source("code/functions/prepare_env.R")

res_dir <- "results/networks/"
fig_dir <- check_dir("results/figures/networks/")

object <- readRDS(
  "../../data/BrainData/processed/GSE97942/GSE97942_cerebellum_processed.rds"
)

p1 <- DimPlot(
  object,
  reduction = "umap",
  cols = colors9,
  group.by = "Sample",
  label = FALSE
) +
  coord_fixed()
p2 <- DimPlot(
  object,
  reduction = "umap.rpca",
  cols = colors9,
  group.by = "Sample",
  label = FALSE
) +
  coord_fixed()
p3 <- p1 + p2 +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")
ggsave(
  file.path(fig_dir, "sample_umap_GSE97942.pdf"),
  p3,
  width = 8, height = 4
)

p4 <- DimPlot(
  object,
  reduction = "umap",
  cols = color_celltype,
  group.by = "CellType",
  pt.size = 0.5,
  label = FALSE
) +
  coord_fixed()
p5 <- DimPlot(
  object,
  reduction = "umap.rpca",
  cols = color_celltype,
  group.by = "CellType",
  pt.size = 0.5,
  label = FALSE
) +
  coord_fixed()
p6 <- p4 + p5 +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")
ggsave(
  file.path(fig_dir, "celltype_umap_GSE97942.pdf"),
  p6,
  width = 6, height = 6
)
