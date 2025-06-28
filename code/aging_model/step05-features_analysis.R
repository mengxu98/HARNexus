library(scop)
library(ClusterGVis)
source("code/functions/aging_plots.R")

fig_dir <- "/Users/mx/Study/repositories/HARNexus/results/aging_model_60_100/models/plots/"

expression_matrix <- read.csv(
  "data/ukb/olink_data_imputed_mean_drop30.csv",
  row.names = 1
)

folders <- c("har", "hagr", "brain")
marker_genes_list <- list(
  har = c(
    "EDA2R",
    "ELN",
    "CDCP1",
    "GFAP",
    "CXCL14",
    "IGDCC4",
    "CDH3",
    "SFRP4",
    "CDON",
    "IL17D"
  ),
  hagr = c(
    "ELN",
    "EFEMP1",
    "FAS",
    "RET",
    "EGFR",
    "ERBB2",
    "SNCG",
    "IGFBP3",
    "GH1",
    "IGF1R"
  ),
  brain = c(
    "NEFL",
    "GFAP",
    "TREM2",
    "BCAN",
    "PTPRZ1",
    "MOG",
    "PTPRR",
    "CA14",
    "C1QL2",
    "MEPE"
  )
)

aggregate_data <- function(
    matrix,
    meta_data,
    num_group = 5,
    add_group_name = TRUE,
    use_age = TRUE) {
  if (use_age) {
    meta_data$Predict_age <- round(meta_data$Predict_age)
    unique_ages <- sort(unique(meta_data$Predict_age))
    matrix_list <- list()
    for (age in unique_ages) {
      group <- meta_data[meta_data$Predict_age == age, ]
      matrix_list[[as.character(age)]] <- matrix[rownames(group), ]
    }
    matrix_list2 <- lapply(matrix_list, function(x) colMeans(x))
    aggregate_matrix <- do.call(rbind, matrix_list2)
    if (add_group_name) {
      rownames(aggregate_matrix) <- paste0("Age_", unique_ages)
    }
  } else {
    age_group <- quantile(meta_data$Predict_age, seq(0, 1, 1 / num_group))
    age_group <- round(age_group, 0)
    matrix_list <- list()
    for (i in 1:(length(age_group) - 1)) {
      if (i == 1) {
        group <- meta_data[meta_data$Predict_age <= age_group[i + 1], ]
      } else {
        group <- meta_data[
          meta_data$Predict_age > age_group[i] &
            meta_data$Predict_age <= age_group[i + 1],
        ]
      }
      matrix_list[[i]] <- matrix[rownames(group), ]
    }
    matrix_list2 <- lapply(matrix_list, function(x) colMeans(x))
    aggregate_matrix <- do.call(rbind, matrix_list2)
    if (add_group_name) {
      rownames(aggregate_matrix) <- paste0("Age_", age_group[-1])
    }
  }
  aggregate_matrix <- na.omit(aggregate_matrix)
  return(aggregate_matrix)
}

for (folder in folders) {
  features_file <- paste0(
    "results/aging_model_60_100/features_selection/",
    folder,
    "/features_optimized.csv"
  )
  age_results_file <- paste0(
    "results/aging_model_60_100/features_selection/",
    folder,
    "/age_results.csv"
  )
  features <- read.csv(features_file)[, 1]
  mark_genes <- marker_genes_list[[folder]]

  meta_data <- read.csv(age_results_file)
  rownames(meta_data) <- meta_data$eid
  meta_data <- meta_data[, c("eid", "Predict_age")]
  meta_data <- na.omit(meta_data)
  meta_data <- meta_data[order(meta_data$Predict_age), ]
  expression_matrix_features <- expression_matrix[, features]
  expression_matrix_sorted <- expression_matrix_features[rownames(meta_data), ]

  meta_data$group <- "UKB"
  obj <- Seurat::CreateSeuratObject(
    t(expression_matrix),
    meta.data = meta_data
  )
  obj <- Seurat::NormalizeData(obj)
  obj <- Seurat::ScaleData(obj)

  ht <- DynamicHeatmap(
    srt = obj,
    lineages = c("Predict_age"),
    features = features,
    features_label = mark_genes,
    use_fitted = TRUE,
    min_expcells = 0,
    r.sq = 0,
    dev.expl = 0,
    padjust = 1,
    n_split = 3,
    nlabel = 60,
    exp_legend_title = "Z-score",
    column_title = "Predict age",
    split_method = "mfuzz",
    label_size = 10,
    label_color = "black",
    db = "GO_BP",
    anno_terms = TRUE,
    heatmap_palette = "viridis",
    feature_annotation_palcolor = list(
      c("forestgreen")
    ),
    pseudotime_label = 25,
    pseudotime_label_color = "gray",
    pseudotime_label_linetype = 1,
    pseudotime_label_linewidth = 1,
    height = 2,
    width = 1.8
  )
  pdf(
    paste0(fig_dir, "Fig.6e_", folder, ".pdf"),
    width = 15,
    height = 5,
    onefile = FALSE
  )
  print(ht$plot)
  dev.off()

  aggregate_matrix <- aggregate_data(
    expression_matrix_sorted,
    meta_data,
    num_group = 20,
    add_group_name = TRUE,
    use_age = TRUE
  )

  expression_matrix_sorted_t <- t(aggregate_matrix)
  expression_matrix_sorted_t <- as.data.frame(scale(expression_matrix_sorted_t))
  getClusters(expression_matrix_sorted_t)
  cm <- clusterData(
    expression_matrix_sorted_t,
    cluster.method = "mfuzz",
    cluster.num = 3
  )
  visCluster(
    object = cm,
    plot.type = "line"
  )
  ck <- clusterData(
    expression_matrix_sorted_t,
    cluster.method = "kmeans",
    cluster.num = 3
  )
  visCluster(
    object = ck,
    plot.type = "line"
  )

  visCluster(
    object = ck,
    plot.type = "both"
  )

  library(org.Hs.eg.db)

  enrich_go_bp <- enrichCluster(
    ck,
    OrgDb = org.Hs.eg.db,
    type = "BP",
    pvalueCutoff = 0.05,
    topn = 5,
    seed = 1
  )
  enrich_kegg <- enrichCluster(
    ck,
    OrgDb = org.Hs.eg.db,
    type = "KEGG",
    pvalueCutoff = 0.05,
    topn = 5,
    seed = 1
  )

  pdf(
    paste0(fig_dir, "protein_heatmap_60_cluster_", folder, ".pdf"),
    width = 12,
    height = 3.5,
    onefile = FALSE
  )
  visCluster(
    object = ck,
    plot.type = "both",
    column_names_rot = 30,
    show_row_dend = FALSE,
    ms.col = c("#0099CC", "grey90", "#CC3333"),
    line.size = 0.01,
    mline.col = "#dd3c23",
    markGenes = mark_genes,
    markGenes.side = "left",
    add.sampleanno = FALSE,
    box.arg = c(0.03, "grey20"),
    line.side = "left",
    term.text.limit = c(8, 10),
    annoTerm.data = enrich_go_bp,
    go.col = rep(jjAnno::useMyCol("calm", n = 3), each = 5),
    genes.gp = c("italic", 11, NA),
    subgroup.anno = c("C1", "C2")
  )
  dev.off()
}


file_name <- "merged_seurat_object_PFC_processed_sub.rds"
if (!file.exists(paste0("results/Brain-scRNA-seq/", file_name))) {
  srt <- readRDS(
    "results/Brain-scRNA-seq/merged_seurat_object_PFC_processed.rds"
  )
  cells <- scop:::select_cells(
    srt, c("S8", "S9"), "developmental_stage"
  )
  srt_sub <- srt[, cells]

  srt_sub <- Seurat::ScaleData(srt_sub)
  srt_sub <- Seurat::RunPCA(srt_sub)
  srt_sub <- Seurat::RunUMAP(srt_sub, dims = 1:20)
  srt_sub <- Seurat::FindNeighbors(srt_sub, dims = 1:20)
  srt_sub <- Seurat::FindClusters(srt_sub, resolution = 0.5)

  srt_sub <- harmony::RunHarmony(
    srt_sub,
    group.by.vars = "dataset",
    reduction = "pca",
    dims.use = 1:30,
    reduction.save = "harmony"
  )

  srt_sub <- Seurat::RunUMAP(
    srt_sub,
    reduction = "harmony",
    dims = 1:30,
    reduction.name = "umap.harmony"
  )

  srt_sub <- Seurat::FindNeighbors(
    srt_sub,
    reduction = "harmony",
    dims = 1:30
  )

  saveRDS(
    srt_sub,
    paste0("results/Brain-scRNA-seq/", file_name)
  )
}

srt_sub <- readRDS(
  paste0("results/Brain-scRNA-seq/", file_name)
)

p0 <- scop::CellDimPlot(
  srt = srt_sub,
  group.by = "cell_name",
  reduction = "umap.harmony",
  theme_use = "theme_bw"
)
ggplot2::ggsave(
  paste0(fig_dir, "p0.pdf"),
  p0,
  width = 4,
  height = 3
)

p1 <- scop::FeatureDimPlot(
  srt = srt_sub,
  features = c("EDA2R", "ELN", "CDCP1", "GFAP", "IL17D", "CXCL14"),
  reduction = "umap.harmony",
  theme_use = "theme_bw"
)
ggplot2::ggsave(
  paste0(fig_dir, "p1.pdf"),
  p1,
  width = 7,
  height = 5
)

ht_60 <- scop::GroupHeatmap(
  srt = srt_sub,
  features = mark_genes,
  group.by = c("developmental_stage", "cell_name"),
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_side = "bottom",
  exp_legend_title = "Z-score",
  column_title = "Cell type",
  column_names_rot = 30,
  width = 2,
  height = 2.5
)
pdf(
  paste0(fig_dir, "p-heatmap_har_60.pdf"),
  width = 8,
  height = 7
)
print(ht_60$plot)
dev.off()

ht <- scop::GroupHeatmap(
  srt = srt_sub,
  features = c(
    "EDA2R", "ELN", "CDCP1", "GFAP", "IL17D", "CXCL14"
  ),
  feature_split = NULL,
  feature_split_by = NULL,
  group.by = "cell_name",
  show_row_names = TRUE,
  show_column_names = TRUE,
  flip = TRUE,
  column_names_side = "bottom",
  exp_legend_title = "Z-score",
  row_title = "Cell type",
  row_title_rot = 90,
  column_names_rot = 30
)
pdf(
  paste0(fig_dir, "Fig.6f.pdf"),
  width = 3.8,
  height = 4
)
print(ht$plot)
dev.off()


ht_10 <- scop::GroupHeatmap(
  srt = srt_sub,
  features = c(
    "EDA2R", "ELN", "CDCP1", "GFAP", "IL17D", "CXCL14", # to high
    "IGDCC4", "CDH3", "SFRP4", "CDON" # to low
  ),
  feature_split = c(rep("C3", 6), rep("C1", 4)),
  feature_split_by = NULL,
  group.by = "cell_name",
  show_row_names = TRUE,
  show_column_names = TRUE,
  flip = TRUE,
  column_names_side = "bottom",
  exp_legend_title = "Z-score",
  row_title = "Cell type",
  row_title_rot = 90,
  column_names_rot = 30,
  height = 2.5,
  width = 3
)
pdf(
  paste0(fig_dir, "Fig.6f_10.pdf"),
  width = 7,
  height = 7
)
print(ht_10$plot)
dev.off()


ht_10_90 <- scop::GroupHeatmap(
  srt = srt_sub,
  features = c(
    "EDA2R", "ELN", "CDCP1", "GFAP", "IL17D", "CXCL14", # to high
    "IGDCC4", "CDH3", "SFRP4", "CDON" # to low
  ),
  feature_split = c(rep("C3", 6), rep("C1", 4)),
  feature_split_by = NULL,
  group.by = "cell_name",
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_side = "bottom",
  exp_legend_title = "Z-score",
  column_title = "Cell type",
  row_title_rot = 90,
  column_names_rot = 45,
  height = 2.5,
  width = 2.2
)
pdf(
  paste0(fig_dir, "Fig.6f_10_90.pdf"),
  width = 7,
  height = 7
)
print(ht_10_90$plot)
dev.off()

p2 <- scop::FeatureStatPlot(
  srt = srt_sub,
  stat.by = c(
    "EDA2R", "ELN", "CDCP1", "GFAP", "IL17D"
  ),
  fill.by = "feature",
  group.by = "cell_name",
  stack = TRUE, flip = TRUE,
  theme_use = "theme_bw"
)
ggplot2::ggsave(
  paste0(fig_dir, "p2.pdf"),
  p2,
  width = 8,
  height = 5
)

p3 <- scop::FeatureStatPlot(
  srt = srt_sub,
  stat.by = c("EDA2R", "ELN", "CDCP1", "GFAP", "IL17D"),
  group.by = "cell_name",
  plot.by = "group",
  fill.by = "feature",
  ncol = 5,
  legend.position = "none",
  theme_use = "theme_bw"
)
ggplot2::ggsave(
  paste0(fig_dir, "p3.pdf"),
  p3,
  width = 12,
  height = 2.5
)
