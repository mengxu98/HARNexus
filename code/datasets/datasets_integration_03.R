rm(list = ls())
gc()

source("code/functions/prepare_env.R")

data_dir <- "../../data/BrainData/integration/"

objects_file <- file.path(data_dir, "objects_processed.rds")
if (!file.exists(objects_file)) {
  objects <- readRDS(file.path(data_dir, "objects.rds"))

  p1 <- VlnPlot(
    objects,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    group.by = "Dataset",
    pt.size = 0, ncol = 3
  )
  ggsave(
    filename = file.path(data_dir, "vlnplot_unintegrated.pdf"),
    plot = p1,
    width = 18,
    height = 5
  )
  objects <- subset(
    objects,
    subset = nFeature_RNA > 1000 &
      nFeature_RNA < 10000 &
      percent.mt < 15
  )

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
  log_message("Objects saved to {.file {objects_file}}")
}
