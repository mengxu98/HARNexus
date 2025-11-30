rm(list = ls())
gc()

source("code/functions/prepare_env.R")

data_dir <- "../../data/BrainData/integration/"

objects_file <- file.path(data_dir, "objects_processed.rds")
if (!file.exists(objects_file)) {
  objects <- readRDS(file.path(data_dir, "objects_raw.rds"))
  remove_patterns <- c(
    "^ERCC", "^RPLP", "^RPSL", "^MT-", "^mt-",
    "^LOC", "^LINC", "^RP[0-9]", "^AC[0-9]", "MALAT1", "^HB[^(P)]",
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
  saveRDS(objects, objects_file)

  object_sub <- subset(objects, subset = Dataset == "GSE97942")
  saveRDS(object_sub, file.path(data_dir, "GSE97942_processed.rds"))

  log_message("Objects saved to {.file {objects_file}}")
}
