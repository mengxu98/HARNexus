library(thisutils)
library(Seurat)

data_dir <- "../../data/BrainData/integration"
res_dir <- "results/networks/har_csn_data/"
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

log_message("Loading objects...")
objects <- readRDS(file.path(data_dir, "objects_celltypes.rds"))

log_message("Processing cerebellum data...")
file_cerebellum <- file.path(data_dir, "objects_cerebellum.rds")
if (!file.exists(file_cerebellum)) {
  objects_cerebellum <- subset(objects, subset = BrainRegion == "Cerebellum")
  saveRDS(objects_cerebellum, file_cerebellum)
}

log_message("Processing metadata...")
meta_data <- objects@meta.data
unique_combinations <- unique(meta_data[, c("BrainRegion", "Stage")])

log_message("Processing data...")
for (i in seq_len(nrow(unique_combinations))) {
  brain_region <- unique_combinations$BrainRegion[i]
  stage <- unique_combinations$Stage[i]
  log_message("Processing {.val {brain_region}} - {.val {stage}}...")
  file_name <- paste0(brain_region, "_", stage, ".rds")
  file_path <- file.path(res_dir, file_name)
  cells_subset <- rownames(meta_data)[meta_data$BrainRegion == brain_region & meta_data$Stage == stage]
  log_message("Found {.val {length(cells_subset)}} cells")
  if (!file.exists(file_path)) {
    objects_subset <- objects[, cells_subset]
    saveRDS(objects_subset, file_path)
  }
  rm(objects_subset)
  gc()
}
