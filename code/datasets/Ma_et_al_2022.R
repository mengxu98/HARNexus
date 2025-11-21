rm(list = ls())
gc()

source("code/functions/prepare_env.R")

data_dir <- "../../data/BrainData/raw/GSE207334"
res_dir <- check_dir("../../data/BrainData/processed/Ma_et_al_2022/")


log_message("Start loading data...")

object2 <- readRDS(
  file.path(data_dir, "Ma_Sestan_mat.rds")
)
object2 <- UpdateSeuratObject(object2)


sample_info <- data.frame(
  Sample = c("HSB106", "HSB189", "HSB340", "HSB628"),
  Age = c(64, 36, 19, 50),
  Sex = c("Male", "Male", "Male", "Female"),
  stringsAsFactors = FALSE
)

metadata <- object2@meta.data
metadata$Cells <- rownames(metadata)
metadata$Dataset <- "Ma_et_al_2022"
metadata$Technology <- "10X Genomics"
metadata$Sequence <- "snRNA-seq"
metadata$Sample <- metadata$samplename
metadata$Sample_ID <- metadata$samplename

metadata <- merge(
  metadata, sample_info,
  by = "Sample", all.x = TRUE
)
rownames(metadata) <- metadata$Cells

mapping <- c(
  "L2-3 IT" = "Excitatory neurons",
  "L3-5 IT-1" = "Excitatory neurons",
  "L3-5 IT-2" = "Excitatory neurons",
  "L3-5 IT-3" = "Excitatory neurons",
  "L5 ET" = "Excitatory neurons",
  "L5-6 NP" = "Excitatory neurons",
  "L6 CT" = "Excitatory neurons",
  "L6 IT-1" = "Excitatory neurons",
  "L6 IT-2" = "Excitatory neurons",
  "L6B" = "Excitatory neurons",
  "KCNG1" = "Excitatory neurons",
  "RB" = "Excitatory neurons",
  "ADARB2" = "Inhibitory neurons",
  "LAMP5 LHX6" = "Inhibitory neurons",
  "LAMP5 RELN" = "Inhibitory neurons",
  "PVALB" = "Inhibitory neurons",
  "PVALB ChC" = "Inhibitory neurons",
  "SST" = "Inhibitory neurons",
  "SST HGF" = "Inhibitory neurons",
  "SST NPY" = "Inhibitory neurons",
  "VIP" = "Inhibitory neurons",
  "Astro" = "Astrocytes",
  "Oligo" = "Oligodendrocytes",
  "OPC" = "Oligodendrocyte progenitor cells",
  "Micro" = "Microglia",
  "Immune" = "Immune cells",
  "Endo" = "Endothelial cells",
  "PC" = "Pericytes",
  "SMC" = "Smooth muscle cells",
  "VLMC" = "Vascular and leptomeningeal cells"
)
metadata$CellType <- mapping[metadata$subclass]
metadata <- metadata[metadata$CellType != "Immune cells", ]

metadata$Brain_Region <- "Dorsolateral prefrontal cortex"
metadata$Region <- "Dorsolateral prefrontal cortex"

column_order <- c(
  "Cells", "Dataset", "Technology", "Sequence", "Sample",
  "Sample_ID", "CellType", "Brain_Region", "Region", "Age", "Sex"
)
metadata <- metadata[, column_order]

# file_snrna <- file.path(res_dir, "snRNA-seq_Human_annot_raw.rds")
# if (!file.exists(file_snrna)) {
#   PrepareEnv()
#   sc <- reticulate::import("scanpy")
#   adata2 <- sc$read_h5ad(
#     file.path(data_dir, "processedData/snRNA-seq_Human_annot.h5ad")
#   )
#   object4 <- adata_to_srt(adata2)
#   saveRDS(object4, file_snrna)
# } else {
#   object4 <- readRDS(file_snrna)
# }

metadata <- na.omit(metadata)

counts <- GetAssayData(object2, layer = "counts")
counts <- counts[, metadata$Cells]
object2 <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata
)

log_message("Save data...")
saveRDS(
  object2,
  file.path(res_dir, "Ma_et_al_2022_processed.rds")
)
