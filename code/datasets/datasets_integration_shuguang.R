res_dir <- "../../data/BrainData/integration/"

library(Seurat)
library(ggplot2)
library(patchwork)
library(future)

options(future.globals.maxSize = Inf)
plan("sequential")

message("Loading objects")
objects <- readRDS(file.path(res_dir, "objects.rds"))

message("Filtering objects")
remove_patterns <- c(
  "^ERCC", "^RPLP", "^RPSL", "^MT-", "^mt-",
  "^LOC", "^LINC", "^RP[0-9]", "^AC[0-9]",
  "^AL[0-9]", "^AP[0-9]", "^CT[0-9]",
  "^CH[0-9]", "^FAM[0-9]", "orf"
)

keep <- !grepl(
  paste(
    remove_patterns,
    collapse = "|"
  ), rownames(objects),
  ignore.case = TRUE
)
objects <- objects[keep, ]

message("Running PCA")
dims <- 1:50

objects <- NormalizeData(objects)
objects <- FindVariableFeatures(objects)
objects <- ScaleData(objects)
objects <- RunPCA(objects)


message("Finding neighbors")
objects <- FindNeighbors(
  objects,
  dims = dims,
  reduction = "pca"
)
message("Finding clusters")
objects <- FindClusters(
  objects,
  resolution = 2,
  cluster.name = "unintegrated_clusters"
)

message("Running UMAP")
objects <- RunUMAP(
  objects,
  dims = dims,
  reduction = "pca",
  reduction.name = "umap.unintegrated"
)

p1 <- DimPlot(
  objects,
  reduction = "umap.unintegrated",
  group.by = c("Dataset", "CellType")
)
ggsave(
  filename = file.path(res_dir, "umap_unintegrated.png"),
  plot = p1,
  width = 10,
  height = 10
)
saveRDS(objects, file.path(res_dir, "objects_unintegrated.rds"))

message("Integrating layers using RPCA")
objects <- IntegrateLayers(
  object = objects,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  verbose = FALSE
)

if (!requireNamespace("harmony", quietly = TRUE)) {
  install.packages("harmony")
}
message("Loading harmony")
library(harmony)
message("Integrating layers using Harmony")
objects <- IntegrateLayers(
  object = objects,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.harmony",
  verbose = FALSE
)

message("Finding neighbors using RPCA")
objects <- FindNeighbors(
  objects,
  reduction = "integrated.rpca", dims = 1:30
)
message("Finding clusters using RPCA")
objects <- FindClusters(
  objects,
  resolution = 2,
  cluster.name = "rpca_clusters"
)

message("Running UMAP using RPCA")
objects <- RunUMAP(
  objects,
  reduction = "integrated.rpca",
  dims = 1:30,
  reduction.name = "umap.rpca"
)

message("Finding neighbors using Harmony")
objects <- FindNeighbors(
  objects,
  reduction = "integrated.harmony", dims = 1:30
)
message("Finding clusters using Harmony")
objects <- FindClusters(
  objects,
  resolution = 2, cluster.name = "harmony_clusters"
)

message("Running UMAP using Harmony")
objects <- RunUMAP(
  objects,
  reduction = "integrated.harmony", dims = 1:30, reduction.name = "umap.harmony"
)

message("Visualizing UMAPs")
p4 <- DimPlot(
  objects,
  reduction = "umap.unintegrated",
  group.by = c("unintegrated_clusters")
)
p5 <- DimPlot(
  objects,
  reduction = "umap.rpca",
  group.by = c("rpca_clusters")
)
p6 <- DimPlot(
  objects,
  reduction = "umap.harmony",
  group.by = c("harmony_clusters")
)
pp <- p4 | p5 | p6
ggsave(
  filename = file.path(res_dir, "umap_integration.png"),
  plot = pp,
  width = 15,
  height = 6
)

message("Visualizing UMAPs by cell type")
p1 <- DimPlot(
  objects,
  reduction = "umap.unintegrated",
  group.by = "CellType"
)
p2 <- DimPlot(
  objects,
  reduction = "umap.rpca",
  group.by = "CellType"
)
p3 <- DimPlot(
  objects,
  reduction = "umap.harmony",
  group.by = "CellType"
)
pp <- p1 | p2 | p3
ggsave(
  filename = file.path(res_dir, "umap_integration_joined.png"),
  plot = pp,
  width = 15,
  height = 6
)

message("Saving LISI data")
pca_raw <- objects@reductions$pca@cell.embeddings
pca_rpca <- objects@reductions$integrated.rpca@cell.embeddings
umap_unintegrated <- objects@reductions$umap.unintegrated@cell.embeddings
umap_rpca <- objects@reductions$umap.rpca@cell.embeddings
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


message("Saving objects")
saveRDS(
  objects,
  file.path(res_dir, "objects_processed.rds")
)

message("Joining layers")
objects <- JoinLayers(objects)
objects <- NormalizeData(objects)
objects <- FindVariableFeatures(objects)
objects <- ScaleData(objects)
objects <- RunPCA(objects)

saveRDS(
  objects,
  file.path(res_dir, "objects_processed_joined.rds")
)

# objects <- FindNeighbors(
#   objects,
#   reduction = "integrated.harmony", dims = 1:30
# )
# objects <- FindClusters(
#   objects,
#   resolution = 2,
#   cluster.name = "harmony_clusters"
# )

# objects <- RunHarmony(
#   objects,
#   group.by.vars = "dataset",
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

# objects <- FindNeighbors(
#   objects,
#   reduction = "integrated.harmony",
#   dims = dims
# )
# objects <- FindNeighbors(
#   objects,
#   reduction = "harmony",
#   dims = dims
# )

# saveRDS(
#   objects,
#   objects_processed_file
# )

# pca_raw <- objects@reductions$pca@cell.embeddings
# pca_harmony <- objects@reductions$harmony@cell.embeddings
# umap_raw <- objects@reductions$umap@cell.embeddings
# umap_harmony <- objects@reductions$umap.harmony@cell.embeddings
# meta_data <- objects@meta.data
# lisi_data <- list(
#   pca_raw,
#   pca_harmony,
#   umap_raw,
#   umap_harmony,
#   meta_data
# )
# saveRDS(
#   lisi_data,
#   file.path(res_dir, "lisi_data.rds")
# )
