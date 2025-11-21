rm(list = ls())
gc()

source("code/functions/prepare_env.R")

res_dir <- check_dir("../../data/BrainData/processed/Li_et_al_2018/")

log_message("Start loading data...")
object <- readRDS(
  file.path(res_dir, "Li_et_al_2018.rds")
)

metadata <- object@meta.data
metadata$Cells <- rownames(metadata)
metadata$Dataset <- "Li_et_al_2018"
metadata$Technology <- "10X Genomics"
metadata$Sequence <- metadata$seq_method
metadata$Sample <- metadata$donor_ID
metadata$Sample_ID <- metadata$orig.ident
mapping <- c(
  "Astro"       = "Astrocytes",
  "Endo"        = "Endothelial cells",
  "ExN"         = "Excitatory neurons",
  "InN"         = "Inhibitory neurons",
  "Microglia"   = "Microglia",
  "Oligo"       = "Oligodendrocytes",
  "OPC"         = "Oligodendrocyte progenitor cells",
  "VSMC"        = "Vascular smooth muscle cells"
)
metadata$CellType <- mapping[metadata$original_name]

metadata$Brain_Region <- metadata$region
metadata$Region <- metadata$subregion
metadata$Age <- metadata$donor_age
metadata$Age <- gsub("yr", "", metadata$Age)
metadata <- metadata[metadata$Age != "Unclassified", ]
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
  file.path(res_dir, "Li_et_al_2018_processed.rds")
)
