rm(list = ls())
gc()

source("code/functions/prepare_env.R")

data_dir <- "../../data/BrainData/raw/GSE207334"
res_dir <- check_dir("../../data/BrainData/processed/GSE207334/")

log_message("Start loading data...")

rna_counts_file <- file.path(data_dir, "GSE207334_Multiome_rna_counts.mtx.gz")
rna_genes_file <- file.path(data_dir, "GSE207334_Multiome_rna_genes.txt.gz")
cell_meta_file <- file.path(data_dir, "GSE207334_Multiome_cell_meta.txt.gz")

rna_counts <- Matrix::readMM(gzfile(rna_counts_file))

rna_genes <- read.table(
  gzfile(rna_genes_file),
  header = FALSE, stringsAsFactors = FALSE
)[, 1]

rownames(rna_counts) <- rna_genes

metadata <- read.table(
  gzfile(cell_meta_file),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1
)

colnames(rna_counts) <- rownames(metadata)

sample_info <- data.frame(
  Sample = c("2RT00374N", "RT00382N", "RT00383N", "RT00385N", "RT00390N"),
  Sample_ID = c("HSB6195", "HSB5871", "HSB8050", "HSB6154", "HSB8073"),
  Age = c(45, 60, 43, 68, 51),
  Sex = c("Male", "Male", "Male", "Female", "Female"),
  stringsAsFactors = FALSE
)

metadata$Cells <- rownames(metadata)
metadata$Dataset <- "GSE207334"
metadata$Technology <- "10X Genomics"
metadata$Sequence <- "snRNA-seq"
metadata$Sample <- metadata$samplename

metadata <- merge(
  metadata, sample_info,
  by = "Sample", all.x = TRUE
)
rownames(metadata) <- metadata$Cells

mapping <- c(
  "L2_3_IT"   = "Excitatory neurons",
  "L3_5_IT_1" = "Excitatory neurons",
  "L3_5_IT_2" = "Excitatory neurons",
  "L3_5_IT_3" = "Excitatory neurons",
  "L5_6_NP"   = "Excitatory neurons",
  "L5_PT"     = "Excitatory neurons",
  "L6_CT"     = "Excitatory neurons",
  "L6_IT_1"   = "Excitatory neurons",
  "L6_IT_2"   = "Excitatory neurons",
  "L6b"       = "Excitatory neurons",
  "ADARB2"      = "Inhibitory neurons",
  "LAMP5_LHX6"  = "Inhibitory neurons",
  "LAMP5_RELN"  = "Inhibitory neurons",
  "PVALB"       = "Inhibitory neurons",
  "PVALB_CHC"   = "Inhibitory neurons",
  "SST"         = "Inhibitory neurons",
  "SST_NPY"     = "Inhibitory neurons",
  "VIP"         = "Inhibitory neurons",
  "TH"          = "Inhibitory neurons",
  "Astro"       = "Astrocytes",
  "Oligo"       = "Oligodendrocytes",
  "OPC"         = "Oligodendrocyte progenitor cells",
  "Micro"       = "Microglia",
  "immune"      = "Immune cells",
  "Endo"        = "Endothelial cells",
  "PC"          = "Pericytes",
  "SMC"         = "Smooth muscle cells",
  "VLMC"        = "Vascular and leptomeningeal cells"
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
metadata <- na.omit(metadata)

rna_counts <- rna_counts[, metadata$Cells]
object1 <- CreateSeuratObject(
  counts = rna_counts,
  meta.data = metadata
)

# file_multiome <- file.path(res_dir, "multiome_annot_raw.rds")
# if (!file.exists(file_multiome)) {
#   PrepareEnv()
#   sc <- reticulate::import("scanpy")
#   adata1 <- sc$read_h5ad(
#     file.path(data_dir, "processedData/multiome_annot.h5ad")
#   )
#   object3 <- adata_to_srt(adata1)
#   saveRDS(object3, file_multiome)
# } else {
#   object3 <- readRDS(file_multiome)
# }

log_message("Save data...")

saveRDS(
  object1,
  file.path(res_dir, "GSE207334_processed.rds")
)
