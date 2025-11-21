options(future.globals.maxSize = 1000 * 1024^2)
library(future)
plan("multisession", workers = 6)

data_dir <- "data/seurat_object_list"
res_dir <- "results/Brain-scRNA-seq"
fig_dir <- "results/plots_paper"

dataset_names <- paste0("D", 1:21)
n_colors_needed <- length(dataset_names)

set.seed(42)
color_palette_cluster_seurat <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
  "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
  "#E5C494", "#B3B3B3", "#A6CEE3", "#1F78B4", "#B2DF8A",
  "#33A02C"
)
names(color_palette_cluster_seurat) <- dataset_names

merged_seurat_file <- paste0(
  res_dir, "/merged_seurat_object_processed_PFC.rds"
)
if (!file.exists(merged_seurat_file)) {
  merged_seurat <- subset(
    merged_seurat,
    subset = Area == "PFC"
  )
  saveRDS(
    merged_seurat,
    merged_seurat_file
  )
} else {
  merged_seurat <- readRDS(
    merged_seurat_file
  )
}

if (!file.exists(paste0(res_dir, "/lisi_results.rds"))) {
  lisi_data <- readRDS(paste0(res_dir, "/lisi_data.rds"))
  before_lisi <- compute_lisi(
    X = lisi_data[[3]],
    meta_data = lisi_data[[6]],
    label_colnames = "dataset"
  )
  after_lisi <- compute_lisi(
    X = lisi_data[[4]],
    meta_data = lisi_data[[6]],
    label_colnames = "dataset"
  )
  lisi_results <- cbind(
    before_lisi,
    after_lisi
  )
  names(lisi_results) <- c("Raw", "Harmony")

  saveRDS(lisi_results, paste0(res_dir, "/lisi_results.rds"))
} else {
  lisi_results <- readRDS(paste0(res_dir, "/lisi_results.rds"))
}

lisi_long <- tidyr::gather(
  lisi_results,
  key = "Method",
  value = "LISI"
)
lisi_long$Method <- factor(lisi_long$Method, levels = c("Raw", "Harmony"))

lisi_boxplot <- ggplot(lisi_long, aes(x = Method, y = LISI)) +
  geom_boxplot(
    aes(fill = Method),
    alpha = 0.7,
    width = 0.5,
    outlier.shape = NA
  ) +
  scale_fill_manual(
    values = c("Raw" = "#E41A1C", "Harmony" = "#377EB8")
  ) +
  labs(
    x = "",
    y = "LISI"
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("Raw", "Harmony")),
    label = "p.signif"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )
ggsave(
  paste0(fig_dir, "/lisi_boxplot.pdf"),
  lisi_boxplot,
  width = 5,
  height = 5
)


metadata <- merged_seurat@meta.data

region_stage_mapping <- data.frame(
  development_stage = unique(metadata$development_stage),
  stringsAsFactors = FALSE
)

region_stage_mapping$brain_region <- stringr::str_extract(
  region_stage_mapping$development_stage, "^[^_]+"
)
region_stage_mapping$dev_stage <- stringr::str_extract(
  region_stage_mapping$development_stage, "S\\d+"
)

print("Region-stage mapping:")
print(region_stage_mapping)

brain_regions <- c(
  "ACC", "CBC", "CTX", "DFC", "EC",
  "FC", "HIP", "IG", "ITC", "MDL",
  "MGE", "MTG", "NCX", "OC", "PC",
  "PFC", "Pons", "SN", "TC", "V1C",
  "VMB"
)
dev_stages <- paste0("S", 1:10)

cell_counts_matrix <- matrix(
  0,
  nrow = length(brain_regions),
  ncol = length(dev_stages),
  dimnames = list(brain_regions, dev_stages)
)

for (i in seq_len(nrow(metadata))) {
  region <- stringr::str_extract(
    metadata$development_stage[i], "^[^_]+"
  )
  stage <- stringr::str_extract(
    metadata$development_stage[i], "S\\d+"
  )
  if (!is.na(region) && !is.na(stage) && region %in% brain_regions && stage %in% dev_stages) {
    cell_counts_matrix[region, stage] <- cell_counts_matrix[region, stage] + 1
  }
}

cell_counts_long <- as.data.frame(cell_counts_matrix) %>%
  tibble::rownames_to_column("Region") %>%
  tidyr::pivot_longer(
    cols = -Region,
    names_to = "Stage",
    values_to = "Count"
  )

cell_counts_long$Region <- factor(
  cell_counts_long$Region,
  levels = brain_regions
)
cell_counts_long$Stage <- factor(
  cell_counts_long$Stage,
  levels = dev_stages
)

stage_labels <- c(
  "S1",
  "S2",
  "S3",
  "S4",
  "S5",
  "S6",
  "S7",
  "S8",
  "S9",
  "S10"
)
names(stage_labels) <- dev_stages

stage_descriptions <- c(
  "embryonic (S1):  4 <= age < 8 PCW",
  "early fetal (S2): 8 <= age < 13 PCW",
  "early mid fetal (S3): 13 <= age < 19 PCW",
  "late mid fetal (S4): 19 <= age < 24 PCW",
  "late fetal (S5): 24 <= age < 38 PCW",
  "childhood (S6): 0 <= age < 12 years",
  "adolescence (S7): 12 <= age < 20 years",
  "young adulthood (S8): 20 <= age < 40 years",
  "middle adulthood (S9): 40 <= age < 60 years",
  "late adulthood (S10):  age >= 60 years"
)

heatmap_plot <- ggplot(
  cell_counts_long,
  aes(x = Stage, y = Region, fill = Count)
) +
  labs(
    x = "Development stage",
    y = "Brain region",
    fill = "Cell count"
  ) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradientn(
    colors = c(
      "#FFFFFF",
      "#E6E6FF",
      "#B3C7FF",
      "#6B8FFF",
      "#3366FF",
      "#0033CC",
      "#001F7A"
    ),
    values = c(0, 0.001, 0.1, 0.3, 0.5, 0.7, 1),
    breaks = c(0, 5000, 20000, 50000),
    labels = c("0", "5k", "20k", "50k"),
    name = "Cell count",
    limits = c(0, max(cell_counts_long$Count)),
    na.value = "white"
  ) +
  scale_x_discrete(labels = stage_labels) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.key.width = unit(0.3, "cm"),
    legend.key.height = unit(1, "cm")
  ) +
  coord_fixed()

plot_dataset_before <- DimPlot(
  merged_seurat,
  cols = color_palette_cluster_seurat,
  label = FALSE,
  reduction = "umap",
  group.by = "dataset",
  order = TRUE
) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  ggtitle("Before integration") +
  theme_bw() +
  coord_fixed()

plot_dataset_after <- DimPlot(
  merged_seurat,
  cols = color_palette_cluster_seurat,
  label = FALSE,
  reduction = "umap.harmony",
  group.by = "dataset",
  order = TRUE
) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  ggtitle("Harmony integration") +
  theme_bw() +
  coord_fixed() +
  theme(
    legend.position = "none"
  )

marker_genes <- c(
  "GFAP", "CLDN5", "NRGN",
  "GAD1", "SOX2", "MBP",
  "PDGFRA", "PDGFRB", "P2RY12"
)

color_palette_cluster <- c(
  "Astro" = "#005ea3",
  "Endo" = "#24B700",
  "Micro" = "#00C1AB",
  "OPC" = "#00ACFC",
  "ExN" = "#f0749d",
  "InN" = "#c21b90",
  "NPC" = "#e29828",
  "Olig" = "#5865d3",
  "Perc" = "#c08f09"
)

cell_counts <- table(merged_seurat$cell_name)
merged_seurat$cluster_num <- paste0(
  merged_seurat$cell_name,
  "\n(",
  cell_counts[match(merged_seurat$cell_name, names(cell_counts))],
  ")"
)
color_palette_cluster <- setNames(
  color_palette_cluster,
  paste0(
    names(color_palette_cluster),
    "\n(",
    cell_counts[names(color_palette_cluster)],
    ")"
  )
)

cluster_plots <- DimPlot(
  merged_seurat,
  cols = color_palette_cluster,
  reduction = "umap.harmony",
  group.by = "cluster_num",
  label = TRUE,
  label.size = 3
) +
  labs(x = "UMAP_1", y = "UMAP_2", title = "Cell type") +
  theme_bw() +
  theme(legend.position = "none") +
  coord_fixed()

feature_plots_after <- FeaturePlot(
  merged_seurat,
  features = marker_genes,
  ncol = 9,
  reduction = "umap.harmony",
  cols = c("lightgrey", "#0d0db3cb"),
  order = TRUE
) &
  theme_bw() &
  labs(x = "UMAP_1", y = "UMAP_2") &
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) &
  coord_fixed()

heatmap_plot_dataset <- heatmap_plot +
  plot_dataset_before +
  plot_dataset_after +
  cluster_plots +
  plot_layout(
    ncol = 4
  )
feature_plots_after <- feature_plots_after +
  theme(
    legend.position = "right"
  )

cluster_plots_feature <- lisi_boxplot +
  feature_plots_after +
  plot_layout(widths = c(0.05, 0.95))

all_plots <- heatmap_plot_dataset / cluster_plots_feature +
  plot_layout(heights = c(0.75, 0.25)) +
  plot_annotation(tag_levels = "a")

ggsave(
  paste0(fig_dir, "/fig.S5.pdf"),
  all_plots,
  width = 18,
  height = 7
)

ggsave(
  paste0(fig_dir, "/fig.S5.png"),
  all_plots,
  width = 15,
  height = 6,
  dpi = 600,
  bg = "white"
)


if (FALSE) {
  stage_desc_plot <- ggplot() +
    annotate("text",
      x = 0, y = 0.5,
      label = paste(stage_descriptions, collapse = "\n"),
      hjust = 0.5, vjust = 0.5, size = 3.5
    ) +
    theme_void() +
    theme(
      plot.margin = margin(0, 0, 0, 0),
      plot.background = element_rect(fill = "white", color = NA)
    )

  ggsave(
    paste0(fig_dir, "/stage_desc_plot.pdf"),
    stage_desc_plot,
    width = 7,
    height = 5
  )
}
