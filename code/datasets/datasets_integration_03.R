
source("code/functions/prepare_env.R")

data_dir <- "../../data/BrainData/integration/"

objects_file <- file.path(data_dir, "objects.rds")
if (!file.exists(objects_file)) {
  objects <- readRDS(file.path(data_dir, "objects.rds"))

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
    group.by.vars = c("Dataset", "Sample"),
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
    objects_file
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
