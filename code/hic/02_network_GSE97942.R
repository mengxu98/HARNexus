source("code/functions/prepare_env.R")
source("code/functions/utils_network.R")

res_dir <- check_dir("results/gse97942/")

tfs_raw <- read.csv("results/har_tf/tfs.csv")[, 1]

object <- readRDS(
  "../../data/BrainData/processed/GSE97942/GSE97942_cerebellum_processed.rds"
)
gene_sets_list <- get_celltype_specific_genes(object)

expr_data <- GetAssayData(object, assay = "RNA", layer = "data")

cell_type_vec <- object$CellType

file_path <- file.path(res_dir, "astro_4methods_networks.rds")
if (!file.exists(file_path)) {
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
  write.csv(
    network_genie3,
    file.path(res_dir, "genie3_network.csv"),
    row.names = FALSE
  )
  network_ppcor <- run_ppcor(
    expr,
    regulators = tfs,
    targets = targets
  )
  write.csv(
    network_ppcor,
    file.path(res_dir, "ppcor_network.csv"),
    row.names = FALSE
  )
  network_leap <- run_leap(
    expr,
    regulators = tfs,
    targets = targets
  )
  write.csv(
    network_leap,
    file.path(res_dir, "leap_network.csv"),
    row.names = FALSE
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
  write.csv(
    network_infercsn,
    file.path(res_dir, "infercsn_network.csv"),
    row.names = FALSE
  )
  network_list <- list(
    GENIE3 = network_genie3,
    HARNexus = network_infercsn,
    LEAP = network_leap,
    PPCOR = network_ppcor
  )
  saveRDS(network_list, file_path)
}
