rm(list = ls())
gc()

source("code/functions/prepare_env.R")

data_dir <- "../../data/BrainData/raw/GSE97942/"
res_dir <- check_dir("../../data/BrainData/processed/GSE97942/")

log_message("Start loading data...")
counts_cerebellar_hemisphere <- read.table(
  file.path(data_dir, "CerebellarHem_counts.txt"),
  header = TRUE,
  row.names = 1,
  sep = "\t"
)
counts_frontal_cortex <- read.table(
  file.path(data_dir, "FrontalCortex_counts.txt"),
  header = TRUE,
  row.names = 1,
  sep = "\t"
)

counts_visual_cortex <- read.table(
  file.path(data_dir, "VisualCortex_counts.txt"),
  header = TRUE,
  row.names = 1,
  sep = "\t"
)
common_genes <- Reduce(
  intersect,
  list(
    rownames(counts_cerebellar_hemisphere),
    rownames(counts_frontal_cortex),
    rownames(counts_visual_cortex)
  )
)
counts_cerebellar_hemisphere <- counts_cerebellar_hemisphere[common_genes, ]
counts_frontal_cortex <- counts_frontal_cortex[common_genes, ]
counts_visual_cortex <- counts_visual_cortex[common_genes, ]
counts <- Reduce(
  cbind,
  list(
    counts_cerebellar_hemisphere,
    counts_frontal_cortex,
    counts_visual_cortex
  )
)
rm(
  counts_cerebellar_hemisphere,
  counts_frontal_cortex,
  counts_visual_cortex
)

original_colnames <- colnames(counts)
new_colnames <- sub("^[^_]+_", "", original_colnames)
colnames(counts) <- new_colnames

metadata <- read.csv(
  file.path(data_dir, "metadata.csv"),
  row.names = 1
)
common_cells <- intersect(colnames(counts), rownames(metadata))

counts <- counts[, common_cells, drop = FALSE]
metadata <- metadata[common_cells, , drop = FALSE]

metadata$Cells <- rownames(metadata)
metadata$Dataset <- "GSE97942"
metadata$Technology <- metadata$protocal
metadata$Sequence <- "snRNA-seq"
metadata$Sample <- metadata$orig.ident
metadata$Sample <- gsub("D7_", "", metadata$Sample)
metadata$Sample_ID <- metadata$Sample
mapping <- c(
  "Ast" = "Astrocytes",
  "End" = "Endothelial cells",
  "Ex1" = "Excitatory neurons",
  "Ex2" = "Excitatory neurons",
  "Ex3a" = "Excitatory neurons",
  "Ex3b" = "Excitatory neurons",
  "Ex3c" = "Excitatory neurons",
  "Ex3d" = "Excitatory neurons",
  "Ex3e" = "Excitatory neurons",
  "Ex4" = "Excitatory neurons",
  "Ex5a" = "Excitatory neurons",
  "Ex5b" = "Excitatory neurons",
  "Ex6a" = "Excitatory neurons",
  "Ex6b" = "Excitatory neurons",
  "Ex8" = "Excitatory neurons",
  "Gran" = "Granule cells",
  "In1a" = "Inhibitory neurons",
  "In1b" = "Inhibitory neurons",
  "In1c" = "Inhibitory neurons",
  "In2" = "Inhibitory neurons",
  "In3" = "Inhibitory neurons",
  "In4a" = "Inhibitory neurons",
  "In4b" = "Inhibitory neurons",
  "In6a" = "Inhibitory neurons",
  "In6b" = "Inhibitory neurons",
  "In7" = "Inhibitory neurons",
  "In8" = "Inhibitory neurons",
  "Mic" = "Microglia",
  "Oli" = "Oligodendrocytes",
  "OPC" = "Oligodendrocyte progenitor cells",
  "Per" = "Pericytes",
  "Purk1" = "Purkinje cells",
  "Purk2" = "Purkinje cells"
)
metadata$CellType <- mapping[metadata$priCluster]
metadata$Brain_Region <- metadata$Area
metadata$Region <- metadata$Area
metadata$Age <- gsub(" years", "", metadata$Age)

column_order <- c(
  "Cells", "Dataset", "Technology", "Sequence", "Sample",
  "Sample_ID", "CellType", "Brain_Region", "Region", "Age", "Sex"
)
metadata <- metadata[, column_order]
metadata <- na.omit(metadata)

counts <- counts[, metadata$Cells]
object <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata
)

log_message("Save data...")
saveRDS(
  object,
  file.path(res_dir, "GSE97942_processed.rds")
)
