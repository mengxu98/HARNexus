source("code/functions/prepare_env.R")

data_dir <- "../../data/BrainData/integration/"
fig_dir <- check_dir("figures/datasets")

objects_plot <- readRDS(
  file.path(data_dir, "objects_celltype_plot.rds")
)
objects_plot$Stage <- factor(
  objects_plot$Stage,
  levels = paste0("S", 1:15)
)

Idents(objects_plot) <- "CellType"
celltype_levels <- sort(as.character(unique(objects_plot$CellType)))
objects_plot$CellType <- factor(
  objects_plot$CellType,
  levels = celltype_levels
)

development_stages <- c(
  "Embryonic (4-8 PCW)",
  "Early fetal (8-10 PCW)",
  "Early fetal (10-13 PCW)",
  "Early mid-fetal (13-16 PCW)",
  "Early mid-fetal (16-19 PCW)",
  "Late mid-fetal (19-24 PCW)",
  "Late fetal (24-38 PCW)",
  "Neonatal and early infancy (0-0.5 years)",
  "Late infancy (0.5-1 years)",
  "Early childhood (1-6 years)",
  "Middle and late childhood (6-12 years)",
  "Adolescence (12-20 years)",
  "Young adulthood (20-40 years)",
  "Middle adulthood (40-60 years)",
  "Late adulthood (60+ years)"
)
development_stages2 <- c(
  "Embryonic (4–8 PCW)",
  "Early fetal (8–10 PCW)",
  "Early fetal (10–13 PCW)",
  "Early mid-fetal (13–16 PCW)",
  "Early mid-fetal (16–19 PCW)",
  "Late mid-fetal (19–24 PCW)",
  "Late fetal (24–38 PCW)",
  "Neonatal and early infancy (0–0.5Y)",
  "Late infancy (0.5–1Y)",
  "Early childhood (1–6Y)",
  "Middle and late childhood (6–12Y)",
  "Adolescence (12–20Y)",
  "Young adulthood (20–40Y)",
  "Middle adulthood (40–60Y)",
  "Late adulthood (60+Y)"
)
objects_plot$DevelopmentStage <- factor(
  objects_plot$DevelopmentStage,
  levels = development_stages2
)

stage_def <- data.frame(
  Stage = paste0("S", 1:15),
  DevelopmentStage = development_stages,
  Age_range = c(
    "4-8 PCW", "8-10 PCW", "10-13 PCW", "13-16 PCW", "16-19 PCW",
    "19-24 PCW", "24-38 PCW",
    "0-0.5 years", "0.5-1 years", "1-6 years", "6-12 years",
    "12-20 years", "20-40 years", "40-60 years", "60+ years"
  ),
  stringsAsFactors = FALSE
)
stage_counts <- as.data.frame(
  table(Stage = objects_plot$Stage),
  stringsAsFactors = FALSE
)

names(stage_counts)[2] <- "N_cells"
stage_def <- merge(stage_def, stage_counts, by = "Stage", all.x = TRUE)
stage_def$N_cells[is.na(stage_def$N_cells)] <- 0
stage_def$Stage <- factor(stage_def$Stage, levels = paste0("S", 1:15))
stage_def <- stage_def[order(stage_def$Stage), ]

stage_labels <- as.matrix(
  stage_def[, c("DevelopmentStage", "Age_range", "N_cells")]
)
stage_labels[, "N_cells"] <- format(
  stage_def$N_cells,
  big.mark = ",", trim = TRUE
)
mat_text <- matrix(0, nrow = 15, ncol = 3)
colnames(mat_text) <- c("Development stage", "Age range", "Cell count")
rownames(mat_text) <- as.character(stage_def$Stage)

ra <- rowAnnotation(
  Stage = stage_def$Stage,
  col = list(Stage = color_stages),
  show_legend = FALSE,
  width = unit(1, "mm"),
  annotation_name_rot = 0,
  annotation_name_gp = gpar(fontsize = 9)
)

ht_stage <- Heatmap(
  mat_text,
  left_annotation = ra,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_heatmap_legend = FALSE,
  col = structure("white", names = "0"),
  rect_gp = gpar(col = "gray70", lwd = 0.5),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      stage_labels[i, j],
      x - width / 2, y,
      gp = gpar(fontsize = 9),
      just = "left"
    )
  },
  column_names_side = "bottom",
  column_names_rot = 0,
  column_names_gp = gpar(fontsize = 9),
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 9)
)

pdf(
  file.path(fig_dir, "development_stage_annotation.pdf"),
  width = 8, height = 3
)
draw(ht_stage)
dev.off()


p1 <- CellDimPlot(
  objects_plot,
  reduction = "umap.unintegrated",
  group.by = "Dataset",
  palcolor = colors32,
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)
p2 <- CellDimPlot(
  objects_plot,
  reduction = "umap.rpca",
  group.by = "Dataset",
  palcolor = colors32,
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)

p4 <- DotPlot(
  objects_plot,
  features = marker_genes,
  group.by = "seurat_clusters",
  cols = colors2,
  dot.scale = 5
) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "right",
    legend.key.width = unit(0.2, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.text = element_text(size = 12)
  ) +
  coord_fixed()

ggsave(
  file.path(fig_dir, "dot_plot_markergenes.pdf"),
  p4,
  width = 9,
  height = 20
)

p5 <- FeatureDimPlot(
  objects_plot,
  features = marker_genes,
  reduction = "umap.rpca",
  ncol = 6,
  raster = TRUE,
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)
ggsave(
  file.path(fig_dir, "feature_plots_rpca.pdf"),
  p5,
  width = 13,
  height = 10
)

p6 <- CellDimPlot(
  objects_plot,
  reduction = "umap.rpca",
  group.by = "seurat_clusters",
  palcolor = colors128,
  label = FALSE,
  raster = TRUE,
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)
ggsave(
  file.path(fig_dir, "dim_plots_seurat_clusters_rpca.pdf"),
  p6,
  width = 11, height = 4.5
)

p7 <- CellDimPlot(
  objects_plot,
  reduction = "umap.rpca",
  group.by = "CellType",
  palcolor = color_celltypes[celltype_levels],
  label = FALSE,
  raster = TRUE,
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)

p9 <- CellDimPlot(
  objects_plot,
  reduction = "umap.rpca",
  group.by = "BrainRegion",
  palcolor = colors128,
  label = FALSE,
  raster = TRUE,
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)

color_stages2 <- color_stages
names(color_stages2) <- development_stages2
p10 <- CellDimPlot(
  objects_plot,
  reduction = "umap.rpca",
  group.by = "DevelopmentStage",
  palcolor = color_stages2,
  label = FALSE,
  raster = TRUE,
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)
p11 <- p1 + p9 + p10
ggsave(
  file.path(fig_dir, "brainregion_stage_rpca.pdf"),
  p11,
  width = 28, height = 4
)

p12 <- CellDimPlot(
  objects_plot,
  reduction = "umap.rpca",
  group.by = "Stage",
  palcolor = color_stages,
  label = FALSE,
  raster = TRUE,
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)
p13 <- p2 + p12 + p7
ggsave(
  file.path(fig_dir, "dim_plots_datasets_stage_celltype.pdf"),
  p13,
  width = 19, height = 4
)

marker_genes_split <- data.frame(
  CellType = c(
    rep("Radial glia", 3),
    rep("Endothelial cells", 4),
    rep("Inhibitory neurons", 3),
    rep("Oligodendrocyte progenitor cells", 5),
    rep("Microglia", 3),
    rep("Neuroblasts", 1),
    rep("Excitatory neurons", 3),
    rep("Astrocytes", 5),
    rep("Oligodendrocytes", 3)
  ),
  Genes = c(
    # Radial glia
    "PAX6", "VIM", "GLI3",
    # Endothelial cells
    "CLDN5", "PECAM1", "VWF", "FLT1",
    # Inhibitory neurons
    "GAD1", "GAD2", "SLC6A1",
    # Oligodendrocyte progenitor cells (OPCs)
    "PDGFRA", "CSPG4", "OLIG1", "OLIG2", "SOX10",
    # Microglia
    "CX3CR1", "P2RY12", "CSF1R",
    # Neuroblasts
    "STMN2",
    # Excitatory neurons
    "SLC17A7", "CAMK2A", "SATB2",
    # Astrocytes
    "GFAP", "AQP4", "ALDH1L1", "FGFR3", "GJA1",
    # Oligodendrocytes
    "MOG", "MAG", "CLDN11"
  )
)
celltype_levels <- sort(as.character(unique(marker_genes_split$CellType)))
marker_genes_split$CellType <- factor(
  marker_genes_split$CellType,
  levels = celltype_levels
)
objects_plot$CellType <- factor(
  objects_plot$CellType,
  levels = celltype_levels
)
celltype_colors <- color_celltypes[celltype_levels]
gh <- GroupHeatmap(
  objects_plot,
  exp_legend_title = "Z-score",
  features = marker_genes_split$Genes,
  feature_split = marker_genes_split$CellType,
  group.by = "CellType",
  group_palcolor = celltype_colors,
  cell_annotation_palcolor = celltype_colors,
  feature_split_palcolor = celltype_colors,
  heatmap_palette = "Spectral",
  height = 6,
  width = 3,
  add_dot = TRUE,
  dot_size = unit(6, "mm"),
  nlabel = 0,
  show_row_names = TRUE
)
pdf(
  file.path(fig_dir, "group_heatmap_markergenes.pdf"),
  width = 12,
  height = 8
)
print(gh$plot)
dev.off()
