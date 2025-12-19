rm(list = ls())
gc()

library(inferCSN)
library(Seurat)
library(thisutils)
source("code/functions/celltype_specific_genes.R")

data_dir <- "results/networks/har_csn_data/"
res_dir <- "results/networks/har_csn_atlas/"
res_dir_rds <- file.path(res_dir, "rds/")
res_dir_csv <- file.path(res_dir, "csv/")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(res_dir_rds, showWarnings = FALSE, recursive = TRUE)
dir.create(res_dir_csv, showWarnings = FALSE, recursive = TRUE)

tfs <- read.csv("results/har_tf/tfs.csv")[, 1]

object_files <- list.files(data_dir, pattern = "\\.rds$", full.names = TRUE)

built_csv <- file.path(res_dir, "built_networks_summary.csv")
not_built_csv <- file.path(res_dir, "not_built_networks_summary.csv")

if (file.exists(built_csv)) {
  built_networks <- read.csv(built_csv)
} else {
  built_networks <- data.frame(
    BrainRegion = character(0),
    Stage = character(0),
    CellType = character(0),
    nCells = integer(0),
    nEdges = integer(0),
    nNodes = integer(0),
    nTFs = integer(0),
    nTargets = integer(0)
  )
}

if (file.exists(not_built_csv)) {
  not_built_networks <- read.csv(not_built_csv)
} else {
  not_built_networks <- data.frame(
    BrainRegion = character(0),
    Stage = character(0),
    CellType = character(0),
    Reason = character(0)
  )
}

parse_filename <- function(file_name) {
  file_base <- gsub("\\.rds$", "", file_name)
  stage <- sub(".*_(S[0-9]+)$", "\\1", file_base)
  brain_region <- sub("_(S[0-9]+)$", "", file_base)
  list(BrainRegion = brain_region, Stage = stage)
}

for (file_path in object_files) {
  file_name <- basename(file_path)
  file_base <- gsub("\\.rds$", "", file_name)
  file_info <- parse_filename(file_name)
  brain_region <- file_info$BrainRegion
  stage <- file_info$Stage

  network_file <- file.path(
    res_dir, paste0(file_base, "_network_list.rds")
  )

  if (!file.exists(network_file)) {
    log_message("Reading object: {.val {file_name}}...")
    object <- readRDS(file_path)
    object <- JoinLayers(object)

    expr_mat <- GetAssayData(object, assay = "RNA", layer = "data")

    cell_type_vec <- object$CellType
    cell_types <- sort(unique(cell_type_vec))

    gene_sets_list <- get_celltype_specific_genes(object)

    network_list <- list()
    for (ct in cell_types) {
      log_message("Processing celltype: {.val {ct}}...")
      cells_ct_idx <- which(cell_type_vec == ct)
      n_cells <- length(cells_ct_idx)

      if (n_cells < 50) {
        log_message("Too few cells for {.val {ct}} - {.val {file_name}}")
        not_built_networks <- rbind(
          not_built_networks, data.frame(
            BrainRegion = brain_region,
            Stage = stage,
            CellType = ct,
            Reason = paste0("Too few cells (n=", n_cells, ")")
          )
        )
        next
      }

      expr <- as.matrix(expr_mat[, cells_ct_idx, drop = FALSE])
      expr <- expr[rowSums(expr) > 0, , drop = FALSE]
      targets <- gene_sets_list[[ct]]
      network_list[[ct]] <- tryCatch(
        {
          net <- inferCSN(
            t(expr),
            regulators = tfs,
            targets = targets,
            penalty = "L0L2",
            cross_validation = TRUE,
            n_folds = 5,
            r_squared_threshold = 0.3,
            cores = 200
          )
          built_networks <<- rbind(
            built_networks, data.frame(
              BrainRegion = brain_region,
              Stage = stage,
              CellType = ct,
              nCells = n_cells,
              nEdges = nrow(net),
              nNodes = length(unique(c(net$regulator, net$target))),
              nTFs = length(unique(net$regulator)),
              nTargets = length(unique(net$target))
            )
          )
          write.csv(
            net,
            file.path(res_dir_csv, paste0(file_base, "_", ct, ".csv")),
            row.names = FALSE, quote = FALSE
          )
          saveRDS(
            net,
            file.path(res_dir_rds, paste0(file_base, "_", ct, ".rds"))
          )
          net
        },
        error = function(e) {
          log_message(
            "Error building network for {.val {ct}} - {.val {file_name}}"
          )
          not_built_networks <<- rbind(
            not_built_networks, data.frame(
              BrainRegion = brain_region,
              Stage = stage,
              CellType = ct,
              Reason = paste0("Error: ", conditionMessage(e))
            )
          )
          return(NULL)
        }
      )
    }

    saveRDS(network_list, network_file)
    log_message("Saved network list: {.val {network_file}}")
  }

  write.csv(
    built_networks, built_csv,
    row.names = FALSE, quote = FALSE
  )
  write.csv(
    not_built_networks, not_built_csv,
    row.names = FALSE, quote = FALSE
  )
}
