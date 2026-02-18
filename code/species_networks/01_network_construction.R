source("code/functions/prepare_env.R")
source("code/functions/utils_network.R")


sample_pairs <- list(
  list(human = "h4", chimp = "c4"),
  list(human = "h3", chimp = "c2"),
  list(human = "h1", chimp = "c1")
)

cores <- 6
assay <- "RNA"
layer <- "data"

tfs <- read.csv("results/har_tf/tfs.csv")[, 1]

for (pair in sample_pairs) {
  human_sample <- pair$human
  chimp_sample <- pair$chimp

  log_message(
    "Network construction for {.val {c(human_sample, chimp_sample)}}"
  )

  res_dir <- check_dir(
    paste0("results/species_networks/", human_sample, "_", chimp_sample, "/")
  )
  res_dir_rds <- check_dir(file.path(res_dir, "rds/"))
  res_dir_csv <- check_dir(file.path(res_dir, "csv/"))

  log_message("Processing data...")
  file_human <- file.path(
    res_dir,
    paste0("human_", human_sample, "_object.rds")
  )
  file_chimp <- file.path(
    res_dir,
    paste0("chimp_", chimp_sample, "_object.rds")
  )
  if (!file.exists(file_human) || !file.exists(file_chimp)) {
    objects <- readRDS(
      "../../data/BrainData/processed/GSE192774/GSE192774_processed.rds"
    )

    objects <- objects[, objects$Sample %in% c(human_sample, chimp_sample)]

    object_human <- objects[, objects$Species == "Human" & objects$Sample == human_sample]
    object_human <- NormalizeData(object_human)
    object_human <- FindVariableFeatures(object_human)
    object_human <- ScaleData(object_human)
    object_human <- RunPCA(object_human)
    object_human <- FindNeighbors(object_human)
    object_human <- FindClusters(object_human)
    object_human <- RunUMAP(object_human, dims = 1:20)
    saveRDS(
      object_human,
      file_human
    )

    object_chimp <- objects[, objects$Species == "Chimpanzee" & objects$Sample == chimp_sample]
    object_chimp <- NormalizeData(object_chimp)
    object_chimp <- FindVariableFeatures(object_chimp)
    object_chimp <- ScaleData(object_chimp)
    object_chimp <- RunPCA(object_chimp)
    object_chimp <- FindNeighbors(object_chimp)
    object_chimp <- FindClusters(object_chimp)
    object_chimp <- RunUMAP(object_chimp, dims = 1:20)
    saveRDS(
      object_chimp,
      file_chimp
    )
  } else {
    object_human <- readRDS(file_human)
    object_chimp <- readRDS(file_chimp)
  }

  expr_mat_human <- GetAssayData(
    object_human,
    assay = assay,
    layer = layer
  )
  expr_mat_chimp <- GetAssayData(
    object_chimp,
    assay = assay,
    layer = layer
  )

  log_message("Identifying cell type-specific genes...")
  gene_sets_list_human <- get_celltype_specific_genes(
    object_human,
    assay = assay,
    layer = layer
  )
  gene_sets_list_chimp <- get_celltype_specific_genes(
    object_chimp,
    assay = assay,
    layer = layer
  )
  gene_sets_list <- purrr::map2(
    gene_sets_list_human, gene_sets_list_chimp, function(x, y) {
      intersect(x, y)
    }
  )

  cell_types_human <- sort(unique(object_human$CellType))
  cell_types_chimp <- sort(unique(object_chimp$CellType))

  log_message(
    "Cell types in Human: {.val {cell_types_human}}"
  )
  log_message(
    "Cell types in Chimpanzee: {.val {cell_types_chimp}}"
  )

  human_networks <- list()
  chimp_networks <- list()

  built_networks <- data.frame(
    Species = character(0),
    Sample = character(0),
    CellType = character(0),
    nCells = integer(0),
    nAllGenes = integer(0),
    nCellTypeSpecificGenes = integer(0),
    nEdges = integer(0),
    nNodes = integer(0),
    nTFs = integer(0),
    nTargets = integer(0)
  )

  not_built_networks <- data.frame(
    Species = character(0),
    Sample = character(0),
    CellType = character(0),
    Reason = character(0)
  )

  human_networks_file <- file.path(
    res_dir, paste0("human_", human_sample, "_networks.rds")
  )

  if (file.exists(human_networks_file)) {
    log_message(
      "Reading existing Human networks from {.file {human_networks_file}}..."
    )
    human_networks <- readRDS(human_networks_file)
    for (ct in names(human_networks)) {
      if (!is.null(human_networks[[ct]])) {
        csv_file <- file.path(
          res_dir_csv,
          paste0("human_", human_sample, "_", ct, ".csv")
        )
        if (!file.exists(csv_file)) {
          write.csv(
            human_networks[[ct]],
            csv_file,
            row.names = FALSE,
            quote = FALSE
          )
        }
      }
    }
  } else {
    log_message("Building networks for Human ({.val {human_sample}})...")

    for (ct in cell_types_human) {
      log_message("Processing Human celltype: {.val {ct}}...")
      cells_ct_idx <- which(object_human$CellType == ct)
      n_cells <- length(cells_ct_idx)

      if (n_cells < 50) {
        log_message("Too few cells for Human {.val {ct}} (n={.val {n_cells}})")
        not_built_networks <- rbind(
          not_built_networks,
          data.frame(
            Species = "Human",
            Sample = human_sample,
            CellType = ct,
            Reason = paste0("Too few cells (n=", n_cells, ")")
          )
        )
        next
      }

      expr <- as.matrix(expr_mat_human[, cells_ct_idx, drop = FALSE])
      expr <- expr[rowSums(expr) > 0, , drop = FALSE]
      targets <- gene_sets_list[[ct]]

      tfs_used <- intersect(tfs, rownames(expr))
      targets_used <- intersect(targets, rownames(expr))
      n_all_genes <- length(unique(c(tfs_used, targets_used)))
      n_celltype_specific_genes <- length(targets_used)

      network_file_ct <- file.path(
        res_dir_rds, paste0("human_", human_sample, "_", ct, ".rds")
      )

      if (!file.exists(network_file_ct)) {
        human_networks[[ct]] <- tryCatch(
          {
            net <- inferCSN(
              t(expr),
              regulators = tfs,
              targets = targets,
              penalty = "L0L2",
              cross_validation = TRUE,
              n_folds = 5,
              r_squared_threshold = 0.3,
              cores = cores
            )

            if (is.null(net) || nrow(net) == 0) {
              log_message("  Empty network result for Human {.val {ct}}")
              not_built_networks <<- rbind(
                not_built_networks,
                data.frame(
                  Species = "Human",
                  Sample = human_sample,
                  CellType = ct,
                  Reason = "Empty network result"
                )
              )
              return(NULL)
            }

            built_networks <<- rbind(
              built_networks,
              data.frame(
                Species = "Human",
                Sample = human_sample,
                CellType = ct,
                nCells = n_cells,
                nAllGenes = n_all_genes,
                nCellTypeSpecificGenes = n_celltype_specific_genes,
                nEdges = nrow(net),
                nNodes = length(unique(c(net$regulator, net$target))),
                nTFs = length(unique(net$regulator)),
                nTargets = length(unique(net$target))
              )
            )

            saveRDS(net, file = network_file_ct)
            write.csv(
              net,
              file.path(
                res_dir_csv,
                paste0("human_", human_sample, "_", ct, ".csv")
              ),
              row.names = FALSE,
              quote = FALSE
            )

            net
          },
          error = function(e) {
            log_message(
              "  Error building network for Human {.val {ct}}: {.val {conditionMessage(e)}}"
            )
            not_built_networks <<- rbind(
              not_built_networks,
              data.frame(
                Species = "Human",
                Sample = human_sample,
                CellType = ct,
                Reason = paste0("Error: ", conditionMessage(e))
              )
            )
            return(NULL)
          }
        )
      } else {
        log_message("  Reading existing network from {.file {network_file_ct}}...")
        human_networks[[ct]] <- readRDS(network_file_ct)
      }
    }

    saveRDS(human_networks, human_networks_file)
  }

  chimp_networks_file <- file.path(res_dir, paste0("chimp_", chimp_sample, "_networks.rds"))

  if (file.exists(chimp_networks_file)) {
    log_message(
      "Reading existing Chimpanzee networks from {.file {chimp_networks_file}}..."
    )
    chimp_networks <- readRDS(chimp_networks_file)
    for (ct in names(chimp_networks)) {
      if (!is.null(chimp_networks[[ct]])) {
        csv_file <- file.path(
          res_dir_csv,
          paste0("chimp_", chimp_sample, "_", ct, ".csv")
        )
        if (!file.exists(csv_file)) {
          write.csv(
            chimp_networks[[ct]],
            csv_file,
            row.names = FALSE,
            quote = FALSE
          )
        }
      }
    }
  } else {
    log_message("Building networks for Chimpanzee ({.val {chimp_sample}})...")

    for (ct in cell_types_chimp) {
      log_message("Processing Chimpanzee celltype: {.val {ct}}...")
      cells_ct_idx <- which(object_chimp$CellType == ct)
      n_cells <- length(cells_ct_idx)

      if (n_cells < 50) {
        log_message("Too few cells for Chimpanzee {.val {ct}} (n={.val {n_cells}})")
        not_built_networks <- rbind(
          not_built_networks,
          data.frame(
            Species = "Chimpanzee",
            Sample = chimp_sample,
            CellType = ct,
            Reason = paste0("Too few cells (n=", n_cells, ")")
          )
        )
        next
      }

      expr <- as.matrix(expr_mat_chimp[, cells_ct_idx, drop = FALSE])
      expr <- expr[rowSums(expr) > 0, , drop = FALSE]
      targets <- gene_sets_list[[ct]]

      tfs_used <- intersect(tfs, rownames(expr))
      targets_used <- intersect(targets, rownames(expr))
      n_all_genes <- length(unique(c(tfs_used, targets_used)))
      n_celltype_specific_genes <- length(targets_used)

      network_file_ct <- file.path(
        res_dir_rds, paste0("chimp_", chimp_sample, "_", ct, ".rds")
      )

      if (!file.exists(network_file_ct)) {
        chimp_networks[[ct]] <- tryCatch(
          {
            net <- inferCSN(
              t(expr),
              regulators = tfs,
              targets = targets,
              penalty = "L0L2",
              cross_validation = TRUE,
              n_folds = 5,
              r_squared_threshold = 0.3,
              cores = cores
            )

            if (is.null(net) || nrow(net) == 0) {
              log_message("  Empty network result for Chimpanzee {.val {ct}}")
              not_built_networks <<- rbind(
                not_built_networks,
                data.frame(
                  Species = "Chimpanzee",
                  Sample = chimp_sample,
                  CellType = ct,
                  Reason = "Empty network result"
                )
              )
              return(NULL)
            }

            built_networks <<- rbind(
              built_networks,
              data.frame(
                Species = "Chimpanzee",
                Sample = chimp_sample,
                CellType = ct,
                nCells = n_cells,
                nAllGenes = n_all_genes,
                nCellTypeSpecificGenes = n_celltype_specific_genes,
                nEdges = nrow(net),
                nNodes = length(unique(c(net$regulator, net$target))),
                nTFs = length(unique(net$regulator)),
                nTargets = length(unique(net$target))
              )
            )

            saveRDS(net, file = network_file_ct)
            write.csv(
              net,
              file.path(
                res_dir_csv,
                paste0("chimp_", chimp_sample, "_", ct, ".csv")
              ),
              row.names = FALSE,
              quote = FALSE
            )

            net
          },
          error = function(e) {
            log_message(
              "  Error building network for Chimpanzee {.val {ct}}: {.val {conditionMessage(e)}}"
            )
            not_built_networks <<- rbind(
              not_built_networks,
              data.frame(
                Species = "Chimpanzee",
                Sample = chimp_sample,
                CellType = ct,
                Reason = paste0("Error: ", conditionMessage(e))
              )
            )
            return(NULL)
          }
        )
      } else {
        log_message("  Reading existing network from {.file {network_file_ct}}...")
        chimp_networks[[ct]] <- readRDS(network_file_ct)
      }
    }

    saveRDS(chimp_networks, chimp_networks_file)
  }

  write.csv(
    built_networks,
    file.path(
      res_dir,
      paste0(
        "built_networks_summary_",
        human_sample, "_", chimp_sample, ".csv"
      )
    ),
    row.names = FALSE,
    quote = FALSE
  )

  write.csv(
    not_built_networks,
    file.path(
      res_dir,
      paste0(
        "not_built_networks_summary_",
        human_sample, "_", chimp_sample, ".csv"
      )
    ),
    row.names = FALSE,
    quote = FALSE
  )

  log_message("Network construction completed for {.val {pair}}!")
  log_message("Human networks: {.val {length(human_networks)}} cell types")
  log_message("Chimpanzee networks: {.val {length(chimp_networks)}} cell types")
}
