rm(list = ls())
gc()

source("code/functions/packages.R")

data_dir <- "../../data/BrainData/HARNexus"
res_dir <- "../../data/BrainData/results/HARNexus/"
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

log_message("Start loading data...")

if (!file.exists(paste0(res_dir, "/merged_seurat_object.rds"))) {
  seurat_objects <- list()

  data_files <- list.files(
    data_dir,
    pattern = ".*_seurat\\.Rdata$",
    full.names = TRUE
  )

  all_metadata <- list()
  all_counts <- list()

  for (file in data_files) {
    sample_name <- gsub(".*/(.*?)_seurat\\.Rdata$", "\\1", file)
    load(file)

    if (exists("seurat") && inherits(seurat, "Seurat")) {
      log_message("Processing: {.val {sample_name}}")

      all_counts[[sample_name]] <- GetAssayData(
        seurat,
        layer = "counts"
      )

      meta_data <- seurat@meta.data
      meta_data$sample_id <- sample_name
      all_metadata[[sample_name]] <- meta_data

      rm(seurat)
      gc()
    }
  }

  combined_metadata <- do.call(rbind, all_metadata)

  log_message("Creating merged Seurat object...")
  merged_seurat <- CreateSeuratObject(counts = all_counts)
  merged_seurat@meta.data <- combined_metadata
  rownames(merged_seurat@meta.data) <- Cells(merged_seurat)
  rm(
    all_counts,
    all_metadata,
    combined_metadata
  )
  gc()

  saveRDS(merged_seurat, paste0(res_dir, "/merged_seurat_object.rds"))
} else {
  merged_seurat <- readRDS(paste0(res_dir, "/merged_seurat_object.rds"))
}

if (!file.exists(paste0(res_dir, "/merged_seurat_object_processed.rds"))) {
  merged_seurat <- JoinLayers(merged_seurat)
  merged_seurat <- NormalizeData(merged_seurat)
  merged_seurat <- FindVariableFeatures(merged_seurat)
  merged_seurat <- ScaleData(merged_seurat)
  merged_seurat <- RunPCA(merged_seurat)
  merged_seurat <- FindNeighbors(merged_seurat, dims = 1:50)
  merged_seurat <- FindClusters(merged_seurat, resolution = 2)
  merged_seurat <- RunUMAP(merged_seurat, dims = 1:30)

  merged_seurat <- RunHarmony(
    merged_seurat,
    group.by.vars = "dataset",
    reduction = "pca",
    dims.use = 1:30,
    reduction.save = "harmony"
  )

  merged_seurat <- RunUMAP(
    merged_seurat,
    reduction = "harmony",
    dims = 1:30,
    reduction.name = "umap.harmony"
  )

  merged_seurat <- FindNeighbors(
    merged_seurat,
    reduction = "integrated.harmony",
    dims = 1:30
  )
  merged_seurat <- FindNeighbors(
    merged_seurat,
    reduction = "harmony",
    dims = 1:30
  )

  saveRDS(
    merged_seurat,
    paste0(
      res_dir, "/merged_seurat_object_processed.rds"
    )
  )

  pca_raw <- merged_seurat@reductions$pca@cell.embeddings
  pca_harmony <- merged_seurat@reductions$harmony@cell.embeddings
  umap_raw <- merged_seurat@reductions$umap@cell.embeddings
  umap_harmony <- merged_seurat@reductions$umap.harmony@cell.embeddings
  meta_data <- merged_seurat@meta.data
  lisi_data <- list(
    pca_raw,
    pca_harmony,
    umap_raw,
    umap_harmony,
    meta_data
  )
  saveRDS(
    lisi_data,
    paste0(res_dir, "/lisi_data.rds")
  )
} else {
  merged_seurat <- readRDS(
    paste0(res_dir, "/merged_seurat_object_processed.rds")
  )
}
