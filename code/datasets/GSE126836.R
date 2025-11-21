rm(list = ls())
gc()

source("code/functions/prepare_env.R")

res_dir <- check_dir("../../data/BrainData/processed/GSE126836/")

log_message("Start loading data...")
object <- readRDS(
  file.path(res_dir, "GSE126836.rds")
)

metadata2 <- read.csv(
  file.path(res_dir, "metadata.csv"),
  row.names = 1
)
metadata <- object@meta.data
rownames(metadata) <- gsub("-\\d+$", "", rownames(metadata))
metadata$Cells <- rownames(metadata)
colnames(object) <- gsub("-\\d+$", "", colnames(object))

common_cells <- intersect(
  rownames(metadata),
  rownames(metadata2)
)
metadata <- metadata[common_cells, ]
metadata2 <- metadata2[common_cells, ]

metadata <- cbind(metadata, metadata2)
metadata$Dataset <- "GSE126836"
metadata$Technology <- "10X Genomics"
metadata$Sequence <- "snRNA-seq"
metadata$Sample <- metadata$orig.ident
metadata$Sample_ID <- metadata$sample_ID
cell_name_map <- c(
  Astro = "Astrocytes",
  Endo  = "Endothelial cells",
  ExN   = "Excitatory neurons",
  InN   = "Inhibitory neurons",
  Micro = "Microglia",
  NPC   = "Neural progenitor cells",
  Olig  = "Oligodendrocytes",
  OPC   = "Oligodendrocyte progenitor cells",
  Perc  = "Pericytes"
)
metadata$cell_name_full <- cell_name_map[metadata$cell_name]
metadata$CellType <- metadata$cell_name_full
metadata$Brain_Region <- metadata$region
metadata$Region <- metadata$subregion
metadata$Age <- metadata$donor_age
metadata$Sex <- metadata$donor_gender
metadata$Sex <- ifelse(metadata$Sex == "F", "Female", "Male")

column_order <- c(
  "Cells", "Dataset", "Technology", "Sequence", "Sample",
  "Sample_ID", "CellType", "Brain_Region", "Region", "Age", "Sex"
)
metadata <- metadata[, column_order]
metadata <- na.omit(metadata)

counts <- GetAssayData(object, layer = "counts")
counts <- counts[, metadata$Cells]
object <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata
)

log_message("Save data...")
saveRDS(
  object,
  file.path(res_dir, "GSE126836_processed.rds")
)
