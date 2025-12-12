rm(list = ls())
gc()

source("code/functions/prepare_env.R")

data_dir <- "../../data/BrainData/integration/"
fig_dir <- check_dir("results/figures/datasets")

objects <- readRDS(file.path(data_dir, "objects_processed_integrated.rds"))
p1 <- DimPlot(
  objects,
  reduction = "umap.unintegrated",
  group.by = "Dataset",
  label = TRUE,
  label.size = 8,
  repel = TRUE,
  pt.size = 0.5
)
p2 <- DimPlot(
  objects,
  reduction = "umap.unintegrated",
  group.by = "CellType",
  label = TRUE,
  label.size = 8,
  repel = TRUE,
  pt.size = 0.5
)
p3 <- DimPlot(
  objects,
  reduction = "umap.rpca",
  group.by = "Dataset",
  label = TRUE,
  label.size = 8,
  repel = TRUE,
  pt.size = 0.5
)
p4 <- DimPlot(
  objects,
  reduction = "umap.rpca",
  group.by = "CellType",
  label = TRUE,
  label.size = 8,
  repel = TRUE,
  pt.size = 0.5
)

p <- p1 + p2 + p3 + p4
ggsave(
  file.path(fig_dir, "umap_integrated_cell_type.pdf"),
  p,
  width = 20,
  height = 20
)


metadata <- objects@meta.data
p1 <- StatPlot(
  metadata,
  stat.by = "Dataset",
  group.by = "CellType",
  plot_type = "bar",
  palette = "Chinese"
)
ggsave(
  file.path(fig_dir, "dataset_cell_type_distribution.pdf"),
  p1,
  width = 9,
  height = 5
)

p2 <- StatPlot(
  metadata,
  stat.by = "CellType",
  group.by = "Stage",
  stat_type = "count",
  plot_type = "rose",
  palette = "Chinese",
  theme_use = "theme_blank"
)
ggsave(
  file.path(fig_dir, "cell_type_distribution_stage.pdf"),
  p2,
  width = 9,
  height = 7,
  dpi = 600
)

p3 <- StatPlot(
  metadata,
  stat.by = "Stage",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "rose",
  palette = "Chinese",
  theme_use = "theme_blank"
)
ggsave(
  file.path(fig_dir, "cell_type_distribution_stage.pdf"),
  p3,
  width = 8,
  height = 6
)

p4 <- StatPlot(
  metadata,
  stat.by = "Dataset",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "dot",
  palette = "Chinese"
)
ggsave(
  file.path(fig_dir, "dataset_cell_type_distribution.pdf"),
  p4,
  width = 10,
  height = 10
)

p5 <- StatPlot(
  metadata,
  stat.by = "Dataset",
  group.by = "Stage",
  stat_type = "count",
  plot_type = "dot",
  palette = "Chinese"
)
ggsave(
  file.path(fig_dir, "dataset_stage_cell_type_distribution.pdf"),
  p5,
  width = 10,
  height = 10
)

p6 <- StatPlot(
  metadata,
  stat.by = "Dataset",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "dot",
  palette = "Chinese"
)
ggsave(
  file.path(fig_dir, "dataset_region_cell_type_distribution.png"),
  width = 10,
  height = 10
)

p7 <- StatPlot(
  metadata,
  stat.by = "BrainRegion",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "dot",
  palette = "Chinese"
)
ggsave(
  file.path(fig_dir, "dataset_brain_region_cell_type_distribution.pdf"),
  p7,
  width = 10,
  height = 10
)

p8 <- StatPlot(
  metadata,
  stat.by = c("BrainRegion", "Stage", "CellType"),
  plot_type = "sankey",
  label = TRUE,
  palette = "Chinese"
)
ggsave(
  file.path(fig_dir, "region_stage_cell_type_distribution.pdf"),
  p8,
  width = 18,
  height = 6
)
