rm(list = ls())
gc()

source("code/functions/prepare_env.R")

res_dir <- check_dir("../../data/BrainData/processed/Nowakowski_et_al_2017/")

log_message("Start loading data...")
object <- readRDS(
  file.path(res_dir, "Nowakowski_et_al_2017.rds")
)

metadata <- object@meta.data
metadata$Cells <- rownames(metadata)
metadata$Dataset <- "Nowakowski_et_al_2017"
metadata$Technology <- metadata$seq_tech
metadata$Sequence <- metadata$seq_method
metadata$Sample <- metadata$donor_ID
metadata$Sample_ID <- metadata$orig.ident
mapping <- c(
  "RG-div1" = "Radial glia",
  "RG-div2" = "Radial glia",
  "RG-early" = "Radial glia",
  "oRG" = "Radial glia",
  "tRG" = "Radial glia",
  "vRG" = "Radial glia",
  "MGE-RG1" = "Radial glia",
  "MGE-RG2" = "Radial glia",
  "IPC-div1" = "Intermediate progenitor cells",
  "IPC-div2" = "Intermediate progenitor cells",
  "IPC-nEN1" = "Intermediate progenitor cells",
  "IPC-nEN2" = "Intermediate progenitor cells",
  "IPC-nEN3" = "Intermediate progenitor cells",
  "MGE-IPC1" = "Intermediate progenitor cells",
  "MGE-IPC2" = "Intermediate progenitor cells",
  "MGE-IPC3" = "Intermediate progenitor cells",
  "MGE-div" = "Intermediate progenitor cells",
  "nEN-early1" = "Newborn excitatory neurons",
  "nEN-early2" = "Newborn excitatory neurons",
  "nEN-late" = "Newborn excitatory neurons",
  "EN-PFC1" = "Excitatory neurons",
  "EN-PFC2" = "Excitatory neurons",
  "EN-PFC3" = "Excitatory neurons",
  "EN-V1-1" = "Excitatory neurons",
  "EN-V1-2" = "Excitatory neurons",
  "EN-V1-3" = "Excitatory neurons",
  "nIN1" = "Newborn interneurons",
  "nIN2" = "Newborn interneurons",
  "nIN3" = "Newborn interneurons",
  "nIN4" = "Newborn interneurons",
  "nIN5" = "Newborn interneurons",
  "IN-CTX-CGE1" = "Interneurons",
  "IN-CTX-CGE2" = "Interneurons",
  "IN-CTX-MGE1" = "Interneurons",
  "IN-CTX-MGE2" = "Interneurons",
  "IN-STR" = "Interneurons",
  "Astrocyte" = "Astrocytes",
  "Glyc" = "Astrocytes",
  "OPC" = "Oligodendrocyte progenitor cells",
  "Endothelial" = "Endothelial cells",
  "Mural" = "Vascular mural cells",
  "Microglia" = "Microglia",
  "Choroid" = "Choroid plexus epithelial cells",
  "U1" = "Unclassified cells",
  "U2" = "Unclassified cells",
  "U3" = "Unclassified cells",
  "U4" = "Unclassified cells"
)
metadata$CellType <- mapping[metadata$original_name]
metadata <- metadata[metadata$CellType != "Unclassified cells", ]
metadata$Brain_Region <- metadata$region
metadata$Region <- metadata$region
metadata$Age <- metadata$donor_age
metadata$Age <- gsub("w", " PCW", metadata$Age)
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
  file.path(res_dir, "Nowakowski_et_al_2017_processed.rds")
)
