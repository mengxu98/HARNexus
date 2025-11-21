rm(list = ls())
gc()

library(Seurat)
library(BPCells)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(thisutils)
options(future.globals.maxSize = 3e+09)

res_dir <- "../../data/BrainData/integration/"

objects <- readRDS(file.path(res_dir, "objects_bpcells.rds"))
remove_patterns <- c(
  "^ERCC", "^RPLP", "^RPSL", "^MT-", "^mt-",
  "^LOC", "^LINC", "^RP[0-9]", "^AC[0-9]",
  "^AL[0-9]", "^AP[0-9]", "^CT[0-9]", "^CH[0-9]", "^FAM[0-9]", "orf"
)

keep <- !grepl(
  paste(
    remove_patterns,
    collapse = "|"
  ), rownames(objects),
  ignore.case = TRUE
)
objects <- objects[keep, ]

objects <- FindVariableFeatures(objects)
objects <- SketchData(
  object = objects,
  ncells = 2000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

saveRDS(objects, file.path(res_dir, "objects_bpcells_sketch_2000.rds"))
objects <- readRDS(file.path(res_dir, "objects_bpcells_sketch_2000.rds"))

DefaultAssay(objects) <- "sketch"
objects <- FindVariableFeatures(objects)
objects <- NormalizeData(objects)
gc()
gc()
objects <- ScaleData(objects)
gc()
gc()
objects <- RunPCA(objects)
gc()
gc()
saveRDS(objects, file.path(res_dir, "objects_bpcells_sketch_processed.rds"))
