source("code/functions/prepare_env.R")
source("code/functions/utils_network.R")

res_dir <- check_dir("results/networks/gse97942/")

tfs_raw <- read.csv("results/har_tf/tfs.csv")[, 1]

object <- readRDS(
  "../../data/BrainData/processed/GSE97942/GSE97942_cerebellum_processed.rds"
)
gene_sets_list <- get_celltype_specific_genes(object)

expr_data <- GetAssayData(object, assay = "RNA", layer = "data")

cell_type_vec <- object$CellType
cell_types <- sort(unique(cell_type_vec))

file_path1 <- file.path(res_dir, "celltype_networks.rds")
if (!file.exists(file_path1)) {
  network_celltypes <- list()
  for (cell_type in cell_types) {
    log_message("Processing celltype: {.val {cell_type}}...")
    targets <- gene_sets_list[[cell_type]]
    expr <- expr_data[, which(cell_type_vec == cell_type)]
    expr <- expr[rowSums(expr) > 0, ]
    expr <- Matrix::as.matrix(expr)
    all_genes <- unique(c(targets, tfs_raw))
    all_genes <- intersect(all_genes, rownames(expr))
    expr <- expr[all_genes, ]
    tfs <- intersect(tfs_raw, all_genes)
    targets <- intersect(targets, all_genes)

    network_celltypes[[cell_type]] <- inferCSN::inferCSN(
      Matrix::t(expr),
      regulators = tfs,
      targets = targets,
      penalty = "L0L2",
      cross_validation = TRUE,
      n_folds = 5,
      r_squared_threshold = 0.3,
      cores = 10
    )
  }
  saveRDS(network_celltypes, file_path1)
}

file_path2 <- file.path(res_dir, "astro_4methods_networks.rds")
if (!file.exists(file_path2)) {
  targets <- gene_sets_list[["Astrocytes"]]
  expr <- expr_data[, which(cell_type_vec == "Astrocytes")]
  expr <- expr[rowSums(expr) > 0, ]
  all_genes <- unique(c(targets, tfs_raw))
  all_genes <- intersect(all_genes, rownames(expr))
  expr <- expr[all_genes, ]
  expr <- Matrix::as.matrix(expr)
  tfs <- intersect(tfs_raw, all_genes)
  targets <- intersect(targets, all_genes)

  network_genie3 <- run_genie3(
    expr,
    regulators = tfs,
    targets = targets,
    cores = 10
  )
  network_ppcor <- run_ppcor(
    expr,
    regulators = tfs,
    targets = targets
  )
  network_leap <- run_leap(
    expr,
    regulators = tfs,
    targets = targets
  )
  network_infercsn <- inferCSN::inferCSN(
    Matrix::t(expr),
    regulators = tfs,
    targets = targets,
    penalty = "L0L2",
    cross_validation = TRUE,
    n_folds = 5,
    r_squared_threshold = 0.3,
    cores = 10
  )

  network_list <- list(
    GENIE3 = network_genie3,
    HARNexus = network_infercsn,
    LEAP = network_leap,
    PPCOR = network_ppcor
  )
  saveRDS(network_list, file_path2)
}
