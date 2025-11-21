rm(list = ls())
gc()

source("code/functions/prepare_env.R")

res_dir <- check_dir("../../data/BrainData/processed/GSE103723/")

log_message("Start loading data...")
object <- readRDS(
  file.path(res_dir, "GSE103723.rds")
)

metadata <- object@meta.data
metadata$Cells <- rownames(metadata)
metadata$Dataset <- "GSE103723"
metadata$Technology <- "STRT-seq"
metadata$Sequence <- "scRNA-seq"
metadata$Sample <- metadata$donor_ID
metadata$Sample_ID <- metadata$sample
mapping <- c(
  "Astrocyte" = "Astrocytes",
  "Bergmann glia" = "Bergmann glia",
  "CGE interneuron" = "Interneurons",
  "Cerebellar inhibitory" = "Inhibitory neurons",
  "Committed oligodendrocyte precursor" = "Oligodendrocyte progenitor cells",
  "Deep-layer corticothalamic and 6b" = "Excitatory neurons",
  "Deep-layer intratelencephalic" = "Excitatory neurons",
  "Deep-layer near-projecting" = "Excitatory neurons",
  "Eccentric medium spiny neuron" = "Medium spiny neurons",
  "Ependymal" = "Ependymal cells",
  "Fibroblast" = "Fibroblasts",
  "Hippocampal CA1-3" = "Excitatory neurons",
  "Hippocampal CA4" = "Excitatory neurons",
  "Hippocampal dentate gyrus" = "Excitatory neurons",
  "LAMP5-LHX6 and Chandelier" = "Inhibitory neurons",
  "Lower rhombic lip" = "Rhombic lip cells",
  "MGE interneuron" = "Interneurons",
  "Mammillary body" = "Mammillary body cells",
  "Microglia" = "Microglia",
  "Midbrain-derived inhibitory" = "Inhibitory neurons",
  "Miscellaneous" = "Unclassified cells",
  "Oligodendrocyte" = "Oligodendrocytes",
  "Oligodendrocyte precursor" = "Oligodendrocyte progenitor cells",
  "Splatter" = "Unclassified cells",
  "Thalamic excitatory" = "Excitatory neurons",
  "Upper rhombic lip" = "Rhombic lip cells",
  "Vascular" = "Vascular cells",
  "unannoted" = "Unclassified cells"
)

metadata$CellType <- mapping[metadata$cell_type]
metadata <- metadata[metadata$CellType != "Unclassified cells", ]
metadata$Brain_Region <- metadata$region
metadata$Region <- metadata$region
metadata$Age <- metadata$donor_age
metadata$Age <- ifelse(metadata$Age == "22w", "22 PCW", "23 PCW")
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
  file.path(res_dir, "GSE103723_processed.rds")
)
