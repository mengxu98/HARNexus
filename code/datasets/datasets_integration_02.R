source("code/functions/prepare_env.R")

data_dir <- "../../data/BrainData/integration/"

objects_file <- file.path(data_dir, "objects.rds")
if (!file.exists(objects_file)) {
  objects_list <- readRDS(file.path(data_dir, "objects_list.rds"))
  objects <- merge(
    objects_list[[1]],
    y = objects_list[-1]
  )
  rm(objects_list)
  gc()
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

  saveRDS(objects, objects_file)
}
