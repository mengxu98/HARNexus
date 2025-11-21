rm(list = ls())
gc()

source("code/functions/prepare_env.R")

data_dir <- "../../data/BrainData/raw/HYPOMAP"
res_dir <- check_dir("../../data/BrainData/processed/HYPOMAP/")

log_message("Start loading data...")

object <- readRDS(
  file.path(
    data_dir, "human_HYPOMAP_snRNASeq.rds"
  )
)

metadata <- object@meta.data
metadata <- na.omit(metadata)
metadata$Cells <- rownames(metadata)
metadata$Dataset <- "HYPOMAP"
metadata$Technology <- "10X Genomics"
metadata$Sequence <- "snRNA-seq"
metadata$Sample <- metadata$Donor_ID
metadata$CellType <- metadata$celltype_annotation
metadata$CellType <- sub(".*\\.", "", metadata$CellType)

region_full <- c(
  ARC        = "Arcuate nucleus",
  DMH        = "Dorsomedial hypothalamic nucleus",
  `Fx/OT/ac` = "Fornix region",
  LH         = "Lateral hypothalamus",
  LPOA       = "Lateral preoptic area",
  LTN        = "Lateral tuberal nucleus",
  MAM        = "Mammillary nucleus",
  ME         = "Median eminence",
  MPOA       = "Medial preoptic area",
  Perivent   = "Periventricular nucleus",
  POA        = "Preoptic area",
  PVN        = "Paraventricular nucleus",
  SCN        = "Suprachiasmatic nucleus",
  SON        = "Supraoptic nucleus",
  Thalamaus  = "Thalamus",
  Thalamus   = "Thalamus",
  TMN        = "Tuberomammillary nucleus",
  Vascular   = "Vascular-associated region",
  Vent       = "Ventricular",
  VMH        = "Ventromedial hypothalamic nucleus"
)
metadata$region_full <- region_full[metadata$region]

region_major <- c(
  ARC        = "Hypothalamus",
  DMH        = "Hypothalamus",
  `Fx/OT/ac` = "Hypothalamus",
  LH         = "Hypothalamus",
  LPOA       = "Hypothalamus",
  LTN        = "Hypothalamus",
  MAM        = "Hypothalamus",
  ME         = "Hypothalamus",
  MPOA       = "Hypothalamus",
  Perivent   = "Hypothalamus",
  POA        = "Hypothalamus",
  PVN        = "Hypothalamus",
  SCN        = "Hypothalamus",
  SON        = "Hypothalamus",
  Thalamaus  = "Thalamus",
  Thalamus   = "Thalamus",
  TMN        = "Hypothalamus",
  Vascular   = "Vascular-associated region",
  Vent       = "Ventricular",
  VMH        = "Hypothalamus"
)
metadata$region_major <- region_major[metadata$region]

metadata$Brain_Region <- metadata$region_full
metadata$Region <- metadata$region_major
metadata$Age <- metadata$age_years
metadata$Sex <- metadata$sex

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
  file.path(res_dir, "HYPOMAP_processed.rds")
)
