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

datasets <- c(
  "AllenM1", "EGAD00001006049", "EGAS00001006537",
  "GSE103723", "GSE104276", "GSE144136", "GSE168408", "GSE178175",
  "GSE186538", "GSE199762", "GSE202210", "GSE204683", "GSE207334",
  "GSE212606", "GSE217511", "GSE67835", "GSE81475", "GSE97942",
  "Li_et_al_2018", "Ma_et_al_2022", "Nowakowski_et_al_2017",
  "PRJCA015229", "ROSMAP", "SomaMut"
)

data_dir <- "../../data/BrainData/processed/"
res_dir <- "../../data/BrainData/integration/"

objects_list_file <- file.path(res_dir, "objects_list.rds")
if (!file.exists(objects_list_file)) {
  objects_list <- list()
  for (dataset in datasets) {
    log_message("Loading: {.val {dataset}}")
    file_path <- file.path(
      data_dir, dataset, paste0(dataset, "_processed.rds")
    )
    object <- readRDS(file_path)
    object$unit <- ifelse(grepl("PCW", object$Age), "PCW", "Years")
    object$Stage <- dplyr::case_when(
      object$Age >= 4 & object$Age < 8 & object$unit == "PCW" ~ "Embryonic (4–8 PCW)",
      object$Age >= 8 & object$Age < 10 & object$unit == "PCW" ~ "Early fetal (8–10 PCW)",
      object$Age >= 10 & object$Age < 13 & object$unit == "PCW" ~ "Early fetal (10–13 PCW)",
      object$Age >= 13 & object$Age < 16 & object$unit == "PCW" ~ "Early mid-fetal (13–16 PCW)",
      object$Age >= 16 & object$Age < 19 & object$unit == "PCW" ~ "Early mid-fetal (16–19 PCW)",
      object$Age >= 19 & object$Age < 24 & object$unit == "PCW" ~ "Late mid-fetal (19–24 PCW)",
      object$Age >= 24 & object$Age < 38 & object$unit == "PCW" ~ "Late fetal (24–38 PCW)",
      object$Age >= 0 & object$Age < 0.5 & object$unit == "Years" ~ "Neonatal and early infancy (0–0.5Y)",
      object$Age >= 0.5 & object$Age < 1 & object$unit == "Years" ~ "Late infancy (0.5–1Y)",
      object$Age >= 1 & object$Age < 6 & object$unit == "Years" ~ "Early childhood (1–6Y)",
      object$Age >= 6 & object$Age < 12 & object$unit == "Years" ~ "Middle and late childhood (6–12Y)",
      object$Age >= 12 & object$Age < 20 & object$unit == "Years" ~ "Adolescence (12–20Y)",
      object$Age >= 20 & object$Age < 40 & object$unit == "Years" ~ "Young adulthood (20–40Y)",
      object$Age >= 40 & object$Age < 60 & object$unit == "Years" ~ "Middle adulthood (40–60Y)",
      object$Age >= 60 & object$unit == "Years" ~ "Late adulthood (60+Y)",
      TRUE ~ NA_character_
    )
    object <- JoinLayers(object)
    write_matrix_dir(
      mat = object[["RNA"]]$counts,
      dir = file.path(res_dir, "bpcells", dataset)
    )
    objects_list[[dataset]] <- object
    rm(object)
    gc()
  }
  metadata <- purrr::map_dfr(
    objects_list,
    function(x) {
      meta <- x@meta.data
      for (col in colnames(meta)) {
        meta[[col]] <- as.character(meta[[col]])
      }
      return(meta)
    }
  )
  write.csv(
    metadata,
    file.path(res_dir, "metadata.csv"),
    row.names = FALSE,
    quote = FALSE
  )
  saveRDS(metadata, file.path(res_dir, "metadata.rds"))
  saveRDS(objects_list, objects_list_file)
} else {
  objects_list <- readRDS(objects_list_file)
}

metadata <- readRDS(file.path(res_dir, "metadata.rds"))
counts_list <- list()
for (dataset in datasets) {
  log_message("Loading: {.val {dataset}}")
  counts_mat <- open_matrix_dir(
    dir = file.path(res_dir, "bpcells", dataset)
  )
  counts_list[[dataset]] <- counts_mat
}

objects <- CreateSeuratObject(
  counts = counts_list,
  meta.data = metadata
)
objects
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
saveRDS(
  object = objects,
  file = file.path(res_dir, "objects_bpcells.rds")
)
