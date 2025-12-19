rm(list = ls())
gc()

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

  celltype_map <- c(
    "Ast"   = "Astrocytes",
    "End"   = "Endothelial cells",
    "Gran"  = "Excitatory neurons",
    "Mic"   = "Microglia",
    "Oli"   = "Oligodendrocytes",
    "OPC"   = "Oligodendrocyte progenitor cells",
    "Per"   = "Vascular cells",
    "Purk1" = "Inhibitory neurons",
    "Purk2" = "Inhibitory neurons"
  )

  metadata <- data.frame(
    Cells = colnames(counts),
    CellType = celltype_map[prefix_info],
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
  object <- subset(object, subset = CellType != "Vascular cells")

  object <- split(object, f = object$Sample)
  object <- NormalizeData(object)
  object <- FindVariableFeatures(object)
  object <- ScaleData(object)
  object <- RunPCA(object)
  object <- RunUMAP(object, reduction = "pca", dims = 1:10)
  object <- IntegrateLayers(
    object = object,
    method = RPCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.rpca"
  )
  object <- RunUMAP(
    object,
    reduction = "integrated.rpca",
    dims = 1:10,
    reduction.name = "umap.rpca"
  )
  object <- JoinLayers(object)
  Idents(object) <- "CellType"
  saveRDS(object, file_path)
} else {
  log_message("Cerebellum data already exists: {.val {file_path}}")
}
