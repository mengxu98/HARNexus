rm(list = ls())
gc()

source("code/functions/prepare_env.R")

data_dir <- "../../data/BrainData/integration"
res_dir <- check_dir("results/networks/")
fig_dir <- check_dir("results/figures/networks/")

file_path <- file.path(res_dir, "GSE97942_processed_subset.rds")
if (!file.exists(file_path)) {
  object <- readRDS(file.path(data_dir, "GSE97942_processed.rds"))
  object <- subset(
    object,
    subset = DevelopmentStage == "Adolescence (12–20Y)"
  )
  object <- subset(
    object,
    subset = CellType != "Endothelial cells"
  )
  object <- CreateSeuratObject(
    counts = GetAssayData(object, layer = "counts"),
    meta.data = object@meta.data
  )
  object <- split(object, f = object$orig.ident)
  object <- NormalizeData(object)
  object <- FindVariableFeatures(object)
  object <- ScaleData(object)
  object <- RunPCA(object)
  object <- RunUMAP(object, reduction = "pca", dims = 1:10)
  object <- IntegrateLayers(
    object = object,
    method = RPCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.rpca"
  )
  object <- RunUMAP(
    object,
    reduction = "integrated.rpca",
    dims = 1:10,
    reduction.name = "umap.rpca"
  )
  object <- JoinLayers(object)
  Idents(object) <- "CellType"
  saveRDS(object, file_path)
} else {
  object <- readRDS(file_path)
}

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
  file.path(fig_dir, "01_sample_umap_GSE97942.pdf"),
  p3,
  width = 8, height = 4
)

p4 <- DimPlot(
  object,
  reduction = "umap.rpca",
  cols = color_celltype,
  group.by = "CellType",
  pt.size = 0.5,
  label = FALSE
) +
  coord_fixed()
ggsave(
  file.path(fig_dir, "02_celltype_umap_GSE97942.pdf"),
  p4,
  width = 6, height = 6
)


expr_mat <- GetAssayData(object, assay = "RNA", layer = "data")

cell_type_vec <- object$CellType
cell_types <- sort(unique(cell_type_vec))

min_prop_expressed_base <- 0.10
min_mean_expr <- 0.10
min_sd_expr <- 0.20
min_SI <- 1.20

gene_sets_list <- list()

for (ct in cell_types) {
  log_message("Processing cell type: {.val {ct}}")

  cells_ct_idx <- which(cell_type_vec == ct)
  n_ct <- length(cells_ct_idx)

  if (n_ct < 100) {
    min_prop_expressed <- min(0.05, min_prop_expressed_base)
  } else {
    min_prop_expressed <- min_prop_expressed_base
  }

  expr_ct <- expr_mat[, cells_ct_idx, drop = FALSE]

  prop_expr <- Matrix::rowMeans(expr_ct > 0)
  genes_prop_keep_idx <- which(prop_expr >= min_prop_expressed)

  log_message("  After prop filter: {.val {length(genes_prop_keep_idx)}} genes")
  if (length(genes_prop_keep_idx) == 0) {
    warning("  No genes passed prop filter for ", ct, ", skip.")
    next
  }

  mean_expr_ct_all <- Matrix::rowMeans(expr_ct)
  mean_expr_ct_sel <- mean_expr_ct_all[genes_prop_keep_idx]

  genes_mean_keep_idx <- genes_prop_keep_idx[mean_expr_ct_sel >= min_mean_expr]

  log_message("  After mean filter: {.val {length(genes_mean_keep_idx)}} genes")
  if (length(genes_mean_keep_idx) == 0) {
    warning("  No genes passed mean filter for ", ct, ", skip.")
    next
  }

  expr_ct_dense <- as.matrix(expr_ct[genes_mean_keep_idx, , drop = FALSE])
  sd_expr_ct <- matrixStats::rowSds(expr_ct_dense)

  genes_sd_keep_idx <- genes_mean_keep_idx[sd_expr_ct >= min_sd_expr]

  log_message("  After SD filter: {.val {length(genes_sd_keep_idx)}} genes")
  if (length(genes_sd_keep_idx) == 0) {
    warning("  No genes passed SD filter for ", ct, ", skip.")
    next
  }

  # Cell type specificity index SI
  # SI_g,T = mu_g,T / (max_{others} mu_g,others + epsilon)
  expr_ct_all <- expr_mat[, cells_ct_idx, drop = FALSE]
  mu_ct_all <- Matrix::rowMeans(expr_ct_all)

  cells_other_idx <- which(cell_type_vec != ct)
  expr_other <- expr_mat[, cells_other_idx, drop = FALSE]
  mu_other_all <- Matrix::rowMeans(expr_other)

  epsilon <- 1e-6
  SI_all <- mu_ct_all / (mu_other_all + epsilon)

  SI_sel <- SI_all[genes_sd_keep_idx]
  genes_SI_keep_idx <- genes_sd_keep_idx[SI_sel >= min_SI]

  log_message("  After SI filter: {.val {length(genes_SI_keep_idx)}} genes")
  if (length(genes_SI_keep_idx) == 0) {
    warning("  No genes passed SI filter for ", ct, ", skip.")
    next
  }

  final_genes_ct <- rownames(expr_mat)[genes_SI_keep_idx]
  gene_sets_list[[ct]] <- final_genes_ct

  log_message("  Final genes for {.val {ct}}: {.val {length(final_genes_ct)}}")
}

str(gene_sets_list)

gene_sets_df <- do.call(
  rbind,
  lapply(
    names(gene_sets_list), function(ct) {
      data.frame(
        CellType = ct,
        Gene = gene_sets_list[[ct]],
        stringsAsFactors = FALSE
      )
    }
  )
)

head(gene_sets_df)
table(gene_sets_df$CellType)

write.csv(
  gene_sets_df,
  file = file.path(res_dir, "01_celltype_specific_genes.csv"),
  quote = FALSE,
  row.names = FALSE
)


tf_data <- read.csv("results/step01-har_tf/human/har_tf_pairs.csv")
tfs <- unique(tf_data$TF)
tfs <- intersect(tfs, rownames(expr_mat))

network_list <- list()
for (ct in cell_types) {
  log_message("Processing celltype: {.val {ct}}...")
  cells_ct_idx <- which(cell_type_vec == ct)
  n_ct <- length(cells_ct_idx)
  expr <- expr_mat[, cells_ct_idx]
  expr_dense <- t(as.matrix(expr))
  genes_use <- gene_sets_list[[ct]]
  network_list[[ct]] <- inferCSN(
    expr_dense,
    regulators = tfs,
    targets = genes_use,
    penalty = "L0L2",
    cross_validation = TRUE,
    n_folds = 5,
    r_squared_threshold = 0.3,
    cores = 1
  )
}

lapply(names(network_list), function(x) {
  network <- network_list[[x]]
  log_message(
    "Cell type: {.val {x}}\n",
    "Number of edges: {.val {nrow(network)}}\n",
    "Number of TFs: {.val {length(unique(network$regulator))}}\n",
    "Number of targets: {.val {length(unique(network$target))}}\n"
  )
})

saveRDS(network_list, file.path(res_dir, "02_network_list.rds"))
