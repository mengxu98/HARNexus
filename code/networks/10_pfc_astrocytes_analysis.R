source("code/functions/prepare_env.R")

log_message("Loading data...")
res_dir <- check_dir("results/networks/analysis/")

region_name <- "Prefrontal cortex"
celltype_name <- "Astrocytes"
dims <- 1:10

objects <- readRDS(
  file.path(res_dir, paste0("object_", region_name, "_", celltype_name, ".rds"))
)

selected_stages <- c("S6", "S7", "S11", "S12")

objects <- subset(
  objects,
  subset = Stage %in% selected_stages
)
objects <- subset(objects, subset = Dataset == "GSE168408")

csv_dir <- "results/networks/har_csn_atlas/csv"
network_list <- lapply(
  selected_stages, function(s) {
    f <- file.path(
      csv_dir, paste0("Prefrontal cortex_", s, "_", celltype_name, ".csv")
    )
    if (!file.exists(f)) {
      return(NULL)
    }
    net <- read.csv(f, stringsAsFactors = FALSE)
    regulators <- unique(net$regulator)
    net <- net[!(net$target %in% regulators), ]
    net$Stage <- s
    net
  }
)
network_list <- network_list[!sapply(network_list, is.null)]
network_data <- do.call(rbind, network_list)
target_genes <- unique(network_data$target)
exclude_genes <- c(
  "SCG3", "MEIKIN", "TMEM161B-AS1", "MRVI1",
  "AGT", "DLC1", "CDH20", "SPON1"
)
target_genes <- target_genes[!target_genes %in% exclude_genes]

objects_var <- NormalizeData(objects)
objects_var <- FindVariableFeatures(objects_var, nfeatures = 500)
objects_var <- ScaleData(objects_var)
objects_var <- RunPCA(objects_var)
objects_var <- FindNeighbors(
  objects_var,
  reduction = "pca",
  dims = dims
)
objects_var <- RunUMAP(
  objects_var,
  reduction = "pca",
  dims = dims,
  reduction.name = "umap"
)

objects_var <- FindClusters(
  objects_var,
  resolution = 1
)

objects_target <- NormalizeData(objects)
objects_target <- FindVariableFeatures(objects_target, nfeatures = 500)
objects_target <- ScaleData(objects_target, features = target_genes)
objects_target <- RunPCA(objects_target, features = target_genes)
objects_target <- FindNeighbors(
  objects_target,
  reduction = "pca",
  dims = dims
)
objects_target <- RunUMAP(
  objects_target,
  reduction = "pca",
  dims = dims,
  reduction.name = "umap"
)

objects_target <- FindClusters(
  objects_target,
  resolution = 1
)

cluster2celltype2 <- c(
  "0" = "Late astrocytes",
  "1" = "Late astrocytes",
  "2" = "Late astrocytes",
  "3" = "Late astrocytes",
  "4" = "Late astrocytes",
  "5" = "Middle astrocytes",
  "6" = "Middle astrocytes",
  "7" = "Middle astrocytes",
  "8" = "Early astrocytes",
  "9" = "Early astrocytes",
  "10" = "Early astrocytes",
  "11" = "Middle astrocytes",
  "12" = "Early astrocytes"
)
objects_target$Celltype <- plyr::mapvalues(
  x    = objects_target$seurat_clusters,
  from = names(cluster2celltype2),
  to   = cluster2celltype2
)
Idents(objects_target) <- "Celltype"

objects_var <- objects_var[, colnames(objects_target)]
objects_var$Celltype <- objects_target$Celltype

cluster2celltype1 <- c(
  "0" = "Late astrocytes",
  "1" = "Late astrocytes",
  "2" = "Late astrocytes",
  "3" = "Late astrocytes",
  "4" = "Middle astrocytes",
  "5" = "Middle astrocytes",
  "6" = "Early astrocytes",
  "7" = "Early astrocytes",
  "8" = "Middle astrocytes",
  "9" = "Early astrocytes",
  "10" = "Middle astrocytes",
  "11" = "Early astrocytes",
  "12" = "Late astrocytes"
)
objects_var$Celltype <- plyr::mapvalues(
  x    = objects_var$seurat_clusters,
  from = names(cluster2celltype1),
  to   = cluster2celltype1
)
Idents(objects_var) <- "Celltype"
objects_var$Celltype <- factor(
  objects_var$Celltype,
  levels = c("Early astrocytes", "Middle astrocytes", "Late astrocytes")
)
objects_target$Celltype <- factor(
  objects_target$Celltype,
  levels = c("Early astrocytes", "Middle astrocytes", "Late astrocytes")
)

objects_var <- RunSlingshot(
  objects_var,
  group.by = "Stage",
  reduction = "umap",
  start = "S6"
)
objects_target <- RunSlingshot(
  objects_target,
  group.by = "Stage",
  reduction = "umap",
  start = "S6"
)

objects_var <- RunDynamicFeatures(
  objects_var,
  lineages = "Lineage1",
  features = VariableFeatures(objects_var)
)

objects_target <- RunDynamicFeatures(
  objects_target,
  lineages = "Lineage1",
  features = target_genes
)

saveRDS(
  objects_var,
  file.path(res_dir, "object_pfc_astrocytes_var.rds")
)
saveRDS(
  objects_target,
  file.path(res_dir, "object_pfc_astrocytes_target.rds")
)
