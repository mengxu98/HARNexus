source("code/functions/prepare_env.R")

data_dir <- "../../data/BrainData/raw/GSE97942/"
res_dir <- check_dir("../../data/BrainData/processed/GSE97942/")

file_path <- file.path(res_dir, "GSE97942_cerebellum_processed.rds")
if (!file.exists(file_path)) {
  log_message("Loading cerebellum data...")
  counts <- read.table(
    file.path(data_dir, "CerebellarHem_counts.txt"),
    header = TRUE,
    row.names = 1,
    sep = "\t"
  )

  prefix_info <- sub("_.*", "", colnames(counts))
  print(table(prefix_info))

  colname_parts <- strsplit(colnames(counts), "_")
  sample_info <- sapply(
    colname_parts, function(x) {
      if (length(x) == 4) {
        return(x[3])
      } else if (length(x) >= 3) {
        return(x[2])
      } else if (length(x) >= 2) {
        return(x[2])
      } else {
        return(x[1])
      }
    }
  )

  print(table(sample_info))

  metadata <- data.frame(
    Cells = colnames(counts),
    CellType = NULL,
    CellType_raw = prefix_info,
    Sample = sample_info,
    stringsAsFactors = FALSE,
    row.names = colnames(counts)
  )
  age_info <- c(
    cbm1 = "49",
    cbm2 = "49",
    cbm3 = "48",
    cbm4 = "48",
    cbm5 = "35",
    cbm6 = "35",
    cbm7 = "35",
    cbm8 = "49",
    cbm9 = "49"
  )
  sex_info <- c(
    cbm1 = "Female",
    cbm2 = "Female",
    cbm3 = "Male",
    cbm4 = "Male",
    cbm5 = "Female",
    cbm6 = "Female",
    cbm7 = "Female",
    cbm8 = "Male",
    cbm9 = "Male"
  )
  metadata$Age <- age_info[metadata$Sample]
  metadata$Sex <- sex_info[metadata$Sample]

  object <- CreateSeuratObject(counts, meta.data = metadata)
  object <- JoinLayers(object)
  object <- NormalizeData(object)
  object <- FindVariableFeatures(object)
  object <- ScaleData(object)
  object <- RunPCA(object)
  object <- RunUMAP(object, reduction = "pca", dims = 1:10)
  object <- FindNeighbors(object, reduction = "pca", dims = 1:10)
  object <- FindClusters(object, resolution = 0.5)

  cluster2celltype <- c(
    "0" = "Excitatory neurons",
    "1" = "Excitatory neurons",
    "2" = "Excitatory neurons",
    "3" = "Astrocytes",
    "4" = "Inhibitory neurons",
    "5" = "Inhibitory neurons",
    "6" = "Inhibitory neurons",
    "7" = "Oligodendrocyte progenitor cells",
    "8" = "Oligodendrocytes",
    "9" = "Microglia",
    "10" = "Unclassified",
    "11" = "Endothelial cells",
    "12" = "Unclassified",
    "13" = "Astrocytes"
  )
  object$CellType <- plyr::mapvalues(
    x    = object$seurat_clusters,
    from = names(cluster2celltype),
    to   = cluster2celltype
  )
  Idents(object) <- "CellType"
  object <- subset(object, subset = CellType != "Unclassified")
  object <- subset(object, subset = CellType != "Endothelial cells")

  object <- NormalizeData(object)
  object <- FindVariableFeatures(object)
  object <- ScaleData(object)
  object <- RunPCA(object)
  object <- RunUMAP(object, reduction = "pca", dims = 1:10)
  saveRDS(object, file_path)
} else {
  log_message("Cerebellum data already exists: {.val {file_path}}")
}
