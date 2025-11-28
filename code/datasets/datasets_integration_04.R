rm(list = ls())
gc()

res_dir <- "../../data/BrainData/integration/"

objects_file <- file.path(
  res_dir, "objects_processed_integrated_joined.rds"
)
if (!file.exists(objects_file)) {
  library(Seurat)
  library(patchwork)
  library(future)

  options(future.globals.maxSize = Inf)
  plan(multicore, workers = 64)
  message("Number of parallel workers: ", nbrOfWorkers())

  message("Loading objects")
  objects <- readRDS(file.path(res_dir, "objects_processed.rds"))

  message("Running PCA")
  dims <- 1:50
  n_features <- 3000

  objects <- NormalizeData(objects)
  objects <- FindVariableFeatures(
    objects,
    nfeatures = n_features
  )
  objects <- ScaleData(objects)
  objects <- RunPCA(
    objects,
    features = VariableFeatures(objects)
  )

  message("Finding neighbors")
  objects <- FindNeighbors(
    objects,
    dims = dims,
    reduction = "pca"
  )
  message("Finding clusters")
  objects <- FindClusters(
    objects,
    resolution = 1
  )
  message("Running UMAP")
  objects <- RunUMAP(
    objects,
    dims = dims,
    reduction = "pca",
    reduction.name = "umap.unintegrated"
  )

  message("Integrating layers using RPCA")
  objects <- IntegrateLayers(
    object = objects,
    method = RPCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.rpca"
  )
  message("Running UMAP using RPCA")
  objects <- RunUMAP(
    objects,
    reduction = "integrated.rpca",
    dims = dims,
    reduction.name = "umap.rpca"
  )

  # message("Loading harmony")
  # library(harmony)
  # message("Integrating layers using Harmony")
  # objects <- IntegrateLayers(
  #   object = objects,
  #   method = HarmonyIntegration,
  #   orig.reduction = "pca",
  #   new.reduction = "integrated.harmony"
  # )
  # message("Running UMAP using Harmony")
  # objects <- RunUMAP(
  #   objects,
  #   reduction = "integrated.harmony",
  #   dims = dims,
  #   reduction.name = "umap.harmony"
  # )

  pca_raw <- objects@reductions$pca@cell.embeddings
  pca_rpca <- objects@reductions$integrated.rpca@cell.embeddings
  # pca_harmony <- objects@reductions$integrated.harmony@cell.embeddings
  umap_raw <- objects@reductions$umap.unintegrated@cell.embeddings
  umap_rpca <- objects@reductions$umap.rpca@cell.embeddings
  # umap_harmony <- objects@reductions$umap.harmony@cell.embeddings
  meta_data <- objects@meta.data
  lisi_data <- list(
    pca_raw,
    pca_rpca,
    # pca_harmony,
    umap_raw,
    umap_rpca,
    # umap_harmony,
    meta_data
  )
  saveRDS(
    lisi_data,
    file.path(res_dir, "lisi_data.rds")
  )

  message("Saving objects")
  saveRDS(
    objects,
    file.path(res_dir, "objects_processed_integrated.rds")
  )

  message("Joining layers")
  objects <- JoinLayers(objects)
  objects <- NormalizeData(objects)
  objects <- FindVariableFeatures(
    objects,
    nfeatures = n_features
  )
  objects <- ScaleData(objects)
  objects <- RunPCA(
    objects,
    features = VariableFeatures(objects)
  )
  objects <- FindNeighbors(objects, dims = dims)
  objects <- FindClusters(objects, resolution = 1)
  objects <- RunUMAP(objects, dims = dims)
  # objects <- RunHarmony(
  #   objects,
  #   group.by.vars = c("Dataset", "Sample"),
  #   reduction = "pca",
  #   dims.use = dims,
  #   reduction.save = "harmony"
  # )
  # objects <- RunUMAP(
  #   objects,
  #   reduction = "harmony",
  #   dims = dims,
  #   reduction.name = "umap.harmony"
  # )
  saveRDS(
    objects,
    objects_file
  )
}
