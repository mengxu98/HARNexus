rm(list = ls())
gc()

source("code/functions/prepare_env.R")

data_dir <- "../../data/BrainData/integration/"

objects_file <- file.path(data_dir, "objects.rds")
if (!file.exists(objects_file)) {
  log_message("Loading objects list...")
  objects_list <- readRDS(file.path(data_dir, "objects_list.rds"))
  log_message("Merging objects...")
  objects <- merge(
    objects_list[[1]],
    y = objects_list[-1]
  )
  rm(objects_list)
  gc()

  objects[["percent.mt"]] <- PercentageFeatureSet(
    objects,
    pattern = "^MT-"
  )
  log_message("Saving objects...")
  saveRDS(objects, objects_file)
  log_message("Objects saved to {.file {objects_file}}")
}
