rm(list = ls())
gc()

library(Seurat)
library(BPCells)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(thisutils)
library(future)
library(harmony)

options(future.globals.maxSize = Inf)
plan("sequential")

res_dir <- "../../data/BrainData/integration/"

objects <- readRDS(file.path(res_dir, "objects_bpcells_sketch_processed.rds"))

objects <- IntegrateLayers(
  objects,
  method = HarmonyIntegration,
  orig = "pca",
  new.reduction = "integrated.harmony",
  dims = 1:50,
  k.anchor = 20
)

objects <- IntegrateLayers(
  objects,
  method = RPCAIntegration,
  orig = "pca",
  new.reduction = "integrated.rpca",
  dims = 1:50,
  k.anchor = 20
)

objects <- FindNeighbors(objects, reduction = "integrated.rpca", dims = 1:50)
objects <- FindClusters(objects, resolution = 2)

saveRDS(objects, file.path(res_dir, "objects_bpcells_sketch_integrated.rds"))
