# data: https://brainscope.gersteinlab.org/

rm(list = ls())
gc()

source("code/functions/prepare_env.R")

data_dir <- "../../data/BrainData/raw/BICCN"
res_dir <- check_dir("../../data/BrainData/processed/BICCN/")

log_message("Start loading data...")

counts <- readRDS(
  file.path(data_dir, "BICCN_mat.RDS")
)
metadata <- readRDS(
  file.path(data_dir, "BICCN_meta_share.RDS")
)
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$sample_id

metadata$Dataset_ID <- "BICCN"
metadata$Cells <- rownames(metadata)
metadata$Dataset <- "GSE104276"
metadata$Technology <- "10X Genomics"
metadata$Sequence <- "snRNA-seq"
metadata$Sample <- metadata$donor
metadata$Sample_ID <- metadata$layer
metadata$Brain_Region <- metadata$region
metadata$Region <- "Dorsolateral Prefrontal Cortex"
metadata$Age <- "Unknown"
metadata$Sex <- metadata$sex
metadata$CellType <- metadata$within_area_subclass

mapping <- c(
  ## Excitatory neurons
  "L2/3 IT"      = "Excitatory neurons",
  "L4 IT"        = "Excitatory neurons",
  "L5 IT"        = "Excitatory neurons",
  "L6 IT"        = "Excitatory neurons",
  "L6 IT Car3"   = "Excitatory neurons",
  "L5 ET"        = "Excitatory neurons",
  "L5/6 NP"      = "Excitatory neurons",
  "L6 CT"        = "Excitatory neurons",
  "L6b"          = "Excitatory neurons",
  ## Inhibitory neurons
  "IN"           = "Inhibitory neurons",
  "Vip"          = "Inhibitory neurons",
  "Pvalb"        = "Inhibitory neurons",
  "Sst"          = "Inhibitory neurons",
  "Sst Chodl"    = "Inhibitory neurons",
  "Lamp5"        = "Inhibitory neurons",
  "Lamp5 Lhx6"   = "Inhibitory neurons",
  "Sncg"         = "Inhibitory neurons",
  "Chandelier"   = "Inhibitory neurons",
  "Pax6"         = "Inhibitory neurons",
  ## Glial cells
  "Astro"        = "Astrocytes",
  "Oligo"        = "Oligodendrocytes",
  "OPC"          = "Oligodendrocyte progenitor cells",
  ## Vascular-related cells
  "Endo"         = "Endothelial cells",
  "VLMC"         = "Vascular and leptomeningeal cells",
  "Micro/PVM"    = "Microglia and perivascular macrophages",
  ## Choroid plexus
  "CP"           = "Choroid plexus epithelial cells"
)

metadata$CellType_full <- mapping[metadata$CellType]
metadata$CellType <- metadata$CellType_full

column_order <- c(
  "Cells", "Dataset", "Technology", "Sequence", "Sample",
  "Sample_ID", "CellType", "Brain_Region", "Region", "Age", "Sex"
)

metadata <- metadata[, column_order]
metadata <- na.omit(metadata)

counts <- counts[, metadata$Cells]

object <- Seurat::CreateSeuratObject(
  counts = counts,
  meta.data = metadata
)
saveRDS(
  object,
  file.path(res_dir, "BICCN_processed.rds")
)
