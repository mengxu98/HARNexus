rm(list = ls())
gc()

source("code/functions/prepare_env.R")

data_dir <- "../../data/BrainData/processed/"
res_dir <- check_dir(
  "../../data/BrainData/integration/"
)

log_message("Start loading data...")

datasets <- c(
  "AllenM1", "EGAD00001006049", "EGAS00001006537",
  "GSE103723", "GSE104276", "GSE144136", "GSE168408", "GSE178175",
  "GSE186538", "GSE199762", "GSE202210", "GSE204683", "GSE207334",
  "GSE212606", "GSE217511", "GSE67835", "GSE81475", "GSE97942",
  "Li_et_al_2018", "Ma_et_al_2022", "Nowakowski_et_al_2017",
  "PRJCA015229", "ROSMAP", "SomaMut"
)

objects_processed_file <- file.path(res_dir, "objects_processed.rds")

if (!file.exists(objects_processed_file)) {
  objects_file <- file.path(res_dir, "objects.rds")
  if (!file.exists(objects_file)) {
    # objects_list_file <- file.path(res_dir, "objects_list.rds")
    objects_list_file <- file.path(res_dir, "objects_list_common_genes.rds")
    if (!file.exists(objects_list_file)) {
      objects_list <- list()
      for (dataset in datasets) {
        log_message("Loading: {.val {dataset}}")
        file_path <- file.path(
          data_dir, dataset, paste0(dataset, "_processed.rds")
        )
        objects_list[[dataset]] <- readRDS(file_path)
      }
      saveRDS(objects_list, objects_list_file)
    } else {
      objects_list <- readRDS(objects_list_file)
    }

    for (dataset in datasets) {
      log_message("Processing: {.val {dataset}}")
      print(
        table(objects_list[[dataset]]$CellType)
      )
    }

    for (dataset in datasets) {
      object <- objects_list[[dataset]]
      object$unit <- ifelse(grepl("PCW", object$Age), "PCW", "Years")
      object$Stage <- dplyr::case_when(
        object$Age >= 4 & object$Age < 8 & object$unit == "PCW" ~ "Embryonic (4–8 PCW)",
        object$Age >= 8 & object$Age < 10 & object$unit == "PCW" ~ "Early fetal (8–10 PCW)",
        object$Age >= 10 & object$Age < 13 & object$unit == "PCW" ~ "Early fetal (10–13 PCW)",
        object$Age >= 13 & object$Age < 16 & object$unit == "PCW" ~ "Early mid-fetal (13–16 PCW)",
        object$Age >= 16 & object$Age < 19 & object$unit == "PCW" ~ "Early mid-fetal (16–19 PCW)",
        object$Age >= 19 & object$Age < 24 & object$unit == "PCW" ~ "Late mid-fetal (19–24 PCW)",
        object$Age >= 24 & object$Age < 38 & object$unit == "PCW" ~ "Late fetal (24–38 PCW)",
        object$Age >= 0 & object$Age < 0.5 & object$unit == "Years" ~ "Neonatal & early infancy (0–0.5Y)",
        object$Age >= 0.5 & object$Age < 1 & object$unit == "Years" ~ "Late infancy (0.5–1Y)",
        object$Age >= 1 & object$Age < 6 & object$unit == "Years" ~ "Early childhood (1–6Y)",
        object$Age >= 6 & object$Age < 12 & object$unit == "Years" ~ "Middle & late childhood (6–12Y)",
        object$Age >= 12 & object$Age < 20 & object$unit == "Years" ~ "Adolescence (12–20Y)",
        object$Age >= 20 & object$Age < 40 & object$unit == "Years" ~ "Young adulthood (20–40Y)",
        object$Age >= 40 & object$Age < 60 & object$unit == "Years" ~ "Middle adulthood (40–60Y)",
        object$Age >= 60 & object$unit == "Years" ~ "Late adulthood (60+Y)",
        TRUE ~ NA_character_
      )
      objects_list[[dataset]] <- object
    }

    column_order <- c(
      "Cells", "Dataset", "Technology", "Sequence", "Sample",
      "Sample_ID", "CellType", "Brain_Region", "Region", "Age", "Sex"
    )

    combined_metadata <- purrr::map_dfr(
      objects_list,
      function(x) {
        meta <- x@meta.data[, column_order, drop = FALSE]
        for (col in colnames(meta)) {
          meta[[col]] <- as.character(meta[[col]])
        }
        return(meta)
      }
    )
    # > dim(combined_metadata)
    # [1] 1781695      11

    # for (dataset in datasets) {
    #   object <- objects_list[[dataset]]
    #   adata <- srt_to_adata(object)
    #   adata$write_h5ad(
    #     file.path(res_dir, paste0(dataset, "_processed.h5ad"))
    #   )
    # }


    remove_patterns <- c(
      "^ERCC", "^RPLP", "^RPSL", "^MT-", "^mt-",
      "^LOC", "^LINC", "^RP[0-9]", "^AC[0-9]",
      "^AL[0-9]", "^AP[0-9]", "^CT[0-9]", "^CH[0-9]", "^FAM[0-9]", "orf"
    )

    keep <- !grepl(
      paste(
        remove_patterns,
        collapse = "|"
      ), rownames(objects),
      ignore.case = TRUE
    )
    objects <- objects[keep, ]

    saveRDS(objects, objects_file)
  } else {
    objects <- readRDS(objects_file)
  }

  common_genes <- Reduce(intersect, lapply(objects_list, rownames))
  for (dataset in datasets) {
    objects <- objects_list[[dataset]]
    objects <- objects[common_genes, ]
    objects_list[[dataset]] <- objects
  }
  objects <- merge(
    objects_list[[1]],
    y = objects_list[-1],
    project = "merged"
  )

  remove_patterns <- c(
    "^ERCC", "^RPLP", "^RPSL", "^MT-", "^mt-",
    "^LOC", "^LINC", "^RP[0-9]", "^AC[0-9]",
    "^AL[0-9]", "^AP[0-9]", "^CT[0-9]", "^CH[0-9]", "^FAM[0-9]", "orf"
  )

  keep <- !grepl(
    paste(
      remove_patterns,
      collapse = "|"
    ), rownames(objects),
    ignore.case = TRUE
  )
  objects <- objects[keep, ]
  objects_file2 <- file.path(res_dir, "objects_common_genes.rds")
  saveRDS(objects, objects_file2)

  dims <- 1:50

  objects <- NormalizeData(objects)
  objects <- FindVariableFeatures(objects)
  objects <- ScaleData(objects)
  objects <- RunPCA(objects)

  objects <- FindNeighbors(objects, dims = dims, reduction = "pca")
  objects <- FindClusters(objects, resolution = 2, cluster.name = "unintegrated_clusters")

  obj <- RunUMAP(obj, dims = dims, reduction = "pca", reduction.name = "umap.unintegrated")
  # visualize by batch and cell type annotation
  # cell type annotations were previously added by Azimuth
  DimPlot(obj, reduction = "umap.unintegrated", group.by = c("Method", "predicted.celltype.l2"))


  objects <- RunHarmony(
    objects,
    group.by.vars = "dataset",
    reduction = "pca",
    dims.use = dims,
    reduction.save = "harmony"
  )

  objects <- RunUMAP(
    objects,
    reduction = "harmony",
    dims = dims,
    reduction.name = "umap.harmony"
  )

  objects <- FindNeighbors(
    objects,
    reduction = "integrated.harmony",
    dims = dims
  )
  objects <- FindNeighbors(
    objects,
    reduction = "harmony",
    dims = dims
  )

  saveRDS(
    objects,
    objects_processed_file
  )

  pca_raw <- objects@reductions$pca@cell.embeddings
  pca_harmony <- objects@reductions$harmony@cell.embeddings
  umap_raw <- objects@reductions$umap@cell.embeddings
  umap_harmony <- objects@reductions$umap.harmony@cell.embeddings
  meta_data <- objects@meta.data
  lisi_data <- list(
    pca_raw,
    pca_harmony,
    umap_raw,
    umap_harmony,
    meta_data
  )
  saveRDS(
    lisi_data,
    file.path(res_dir, "lisi_data.rds")
  )
} else {
  objects <- readRDS(
    objects_processed_file
  )
}
