library(Seurat)
library(thisutils)

load_network_data <- function(network_file, celltypes = NULL, verbose = TRUE) {
  log_message(
    "Loading network data from {.val {network_file}}",
    verbose = verbose
  )
  network <- read.csv(network_file, stringsAsFactors = FALSE)

  if (!is.null(celltypes)) {
    network <- network[network$CellType %in% celltypes, ]
    log_message(
      "Filtered network to specified cell types: {.val {length(unique(network$CellType))}} cell types",
      verbose = verbose
    )
  }

  stage_celltype_pairs <- unique(network[, c("Stage", "CellType")])
  log_message(
    "Found {.val {nrow(stage_celltype_pairs)}} unique Stage-CellType pairs",
    verbose = verbose
  )
  return(stage_celltype_pairs)
}

build_region_object <- function(
    brain_region,
    data_dir = "results/networks/har_csn_data/",
    integration_dir = "../../data/BrainData/integration",
    network_file = NULL,
    celltypes = NULL,
    use_cached = TRUE,
    filter_by_network = TRUE,
    verbose = TRUE) {
  stage_celltype_pairs <- NULL
  if (filter_by_network && !is.null(network_file)) {
    if (file.exists(network_file)) {
      stage_celltype_pairs <- load_network_data(
        network_file,
        celltypes = celltypes, verbose = verbose
      )
    } else {
      log_message(
        "Network file not found: {.val {network_file}}, skipping network filtering",
        verbose = verbose
      )
      filter_by_network <- FALSE
    }
  }

  cache_file <- file.path(
    integration_dir, paste0("objects_", brain_region, ".rds")
  )
  cache_file_alt <- file.path(
    integration_dir,
    paste0("objects_", gsub(" ", "_", brain_region), ".rds")
  )

  object <- NULL

  if (use_cached && file.exists(cache_file)) {
    log_message(
      "Loading cached object from {.val {cache_file}}",
      verbose = verbose
    )
    object <- readRDS(cache_file)
  } else if (use_cached && file.exists(cache_file_alt)) {
    log_message(
      "Loading cached object from {.val {cache_file_alt}}",
      verbose = verbose
    )
    object <- readRDS(cache_file_alt)
  }

  if (is.null(object)) {
    region_files <- list.files(
      data_dir,
      pattern = paste0("^", brain_region, "_S[0-9]+\\.rds$"),
      full.names = TRUE
    )

    if (length(region_files) == 0) {
      region_pattern <- gsub(" ", "_", brain_region)
      region_files <- list.files(
        data_dir,
        pattern = paste0("^", region_pattern, "_S[0-9]+\\.rds$"),
        full.names = TRUE
      )
    }

    if (length(region_files) > 0) {
      if (filter_by_network && !is.null(stage_celltype_pairs)) {
        required_stages <- unique(stage_celltype_pairs$Stage)
        region_files <- region_files[sapply(region_files, function(f) {
          stage <- sub(".*_(S[0-9]+)\\.rds$", "\\1", basename(f))
          stage %in% required_stages
        })]
        log_message(
          "Loading {.val {length(region_files)}} stage files (filtered by network)",
          verbose = verbose
        )
      } else {
        log_message(
          "Found {.val {length(region_files)}} stage files for {.val {brain_region}}",
          verbose = verbose
        )
      }

      objects_list <- list()
      for (file_path in region_files) {
        file_name <- basename(file_path)
        stage <- sub(".*_(S[0-9]+)\\.rds$", "\\1", file_name)
        log_message("Loading {.val {stage}}...", verbose = verbose)
        obj_tmp <- readRDS(file_path)

        if (!is.null(celltypes)) {
          meta_tmp <- obj_tmp@meta.data
          keep_cells_tmp <- meta_tmp$CellType %in% celltypes
          if (sum(keep_cells_tmp) > 0) {
            obj_tmp <- obj_tmp[, keep_cells_tmp]
            log_message(
              "  Filtered to {.val {sum(keep_cells_tmp)}} / {.val {length(keep_cells_tmp)}} cells",
              verbose = verbose
            )
            objects_list[[stage]] <- obj_tmp
          } else {
            log_message("  No cells found for specified cell types, skipping", verbose = verbose)
          }
        } else {
          objects_list[[stage]] <- obj_tmp
        }

        rm(obj_tmp)
        gc()
      }

      if (length(objects_list) > 0) {
        log_message(
          "Joining layers for {.val {length(objects_list)}} objects before merging...",
          verbose = verbose
        )
        for (i in seq_along(objects_list)) {
          tryCatch(
            {
              objects_list[[i]] <- JoinLayers(objects_list[[i]])
            },
            error = function(e) {
              log_message(
                "Warning: Failed to join layers for object {.val {i}}: {.val {e$message}}",
                message_type = "warning",
                verbose = verbose
              )
            }
          )
        }

        if (length(objects_list) == 1) {
          object <- objects_list[[1]]
        } else {
          log_message(
            "Merging {.val {length(objects_list)}} objects...",
            verbose = verbose
          )
          object <- merge(
            objects_list[[1]],
            y = objects_list[-1],
            merge.data = TRUE
          )
        }
      }
      object <- JoinLayers(object)
      object <- CreateSeuratObject(
        counts = GetAssayData(object, layer = "counts"),
        meta.data = object@meta.data
      )

      rm(objects_list)
      gc()
    }
  }

  if (filter_by_network && !is.null(stage_celltype_pairs)) {
    log_message("Filtering cells by network data...", verbose = verbose)
    meta <- object@meta.data
    keep_cells <- rep(FALSE, nrow(meta))

    for (i in seq_len(nrow(stage_celltype_pairs))) {
      stage <- stage_celltype_pairs$Stage[i]
      celltype <- stage_celltype_pairs$CellType[i]
      keep_cells <- keep_cells | (meta$Stage == stage & meta$CellType == celltype)
    }

    log_message(
      "Filtering cells: {.val {sum(keep_cells)}} / {.val {ncol(object)}} cells kept",
      verbose = verbose
    )
    object <- object[, keep_cells]
  }

  if (!is.null(object)) {
    object <- JoinLayers(object)
  }
  object <- CreateSeuratObject(
    counts = GetAssayData(object, layer = "counts"),
    meta.data = object@meta.data
  )
  return(object)
}

get_region_summary <- function(object, verbose = TRUE) {
  meta <- object@meta.data

  summary_list <- list(
    cells_per_celltype = table(meta$CellType),
    cells_per_stage_celltype = table(meta$Stage, meta$CellType)
  )

  return(summary_list)
}

preprocess_object <- function(objects) {
  objects <- NormalizeData(objects)
  objects <- FindVariableFeatures(objects)
  objects <- ScaleData(objects)
  objects <- RunPCA(objects)
  objects <- RunUMAP(objects, dims = 1:10)
  return(objects)
}

objects <- build_region_object(
  "Prefrontal cortex",
  celltypes = "Astrocytes",
  network_file = "results/networks/analysis/network_Prefrontal cortex.csv"
)
objects <- preprocess_object(objects)
summary_pfc_astrocytes <- get_region_summary(objects)
print(summary_pfc_astrocytes)
saveRDS(
  objects,
  "results/networks/analysis/object_Prefrontal cortex_Astrocytes.rds"
)

objects <- build_region_object(
  "Prefrontal cortex",
  celltypes = "Excitatory neurons",
  network_file = "results/networks/analysis/network_Prefrontal cortex.csv"
)
objects <- preprocess_object(objects)
summary_pfc_exn <- get_region_summary(objects)
print(summary_pfc_exn)
saveRDS(
  objects,
  "results/networks/analysis/object_Prefrontal cortex_Excitatory neurons.rds"
)

objects <- build_region_object(
  "Cerebral cortex",
  network_file = "results/networks/analysis/network_Cerebral cortex.csv"
)
objects <- preprocess_object(objects)
summary_cc <- get_region_summary(objects)
print(summary_cc)
saveRDS(
  objects,
  "results/networks/analysis/object_Cerebral cortex.rds"
)

objects <- build_region_object(
  "Prefrontal cortex",
  network_file = "results/networks/analysis/network_Prefrontal cortex.csv"
)
objects <- preprocess_object(objects)
summary_pfc <- get_region_summary(objects)
print(summary_pfc)
saveRDS(
  objects,
  "results/networks/analysis/object_Prefrontal cortex.rds"
)

objects <- build_region_object(
  "Prefrontal cortex",
  celltypes = "Astrocytes",
  network_file = "results/networks/analysis/network_Prefrontal cortex.csv"
)
objects <- preprocess_object(objects)
summary_pfc_astrocytes <- get_region_summary(objects)
print(summary_pfc_astrocytes)
saveRDS(
  objects,
  "results/networks/analysis/object_Prefrontal cortex_Astrocytes.rds"
)

objects <- build_region_object(
  "Prefrontal cortex",
  celltypes = "Radial glia",
  network_file = "results/networks/analysis/network_Prefrontal cortex.csv"
)
objects <- preprocess_object(objects)
summary_pfc_radial_glia <- get_region_summary(objects)
print(summary_pfc_radial_glia)
saveRDS(
  objects,
  "results/networks/analysis/object_Prefrontal cortex_Radial glia.rds"
)

objects <- build_region_object(
  "Cerebral cortex",
  celltypes = "Radial glia",
  network_file = "results/networks/analysis/network_Cerebral cortex.csv"
)
objects <- preprocess_object(objects)
summary_cc_radial_glia <- get_region_summary(objects)
print(summary_cc_radial_glia)
saveRDS(
  objects,
  "results/networks/analysis/object_Cerebral cortex_Radial glia.rds"
)
objects <- build_region_object(
  "Cerebral cortex",
  celltypes = "Astrocytes",
  network_file = "results/networks/analysis/network_Cerebral cortex.csv"
)
objects <- preprocess_object(objects)
summary_cc_astrocytes <- get_region_summary(objects)
print(summary_cc_astrocytes)
saveRDS(
  objects,
  "results/networks/analysis/object_Cerebral cortex_Astrocytes.rds"
)
