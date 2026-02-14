source("code/functions/prepare_env.R")

sample_pairs <- list(
  list(human = "h4", chimp = "c4"),
  list(human = "h3", chimp = "c2"),
  list(human = "h1", chimp = "c1")
)

for (pair in sample_pairs) {
  human_sample <- pair$human
  chimp_sample <- pair$chimp

  res_dir <- check_dir(
    paste0("results/species_networks/", human_sample, "_", chimp_sample, "/atac")
  )

  atac_file_sub <- file.path(
    res_dir,
    paste0("GSE192774_ATAC_processed_", human_sample, "_", chimp_sample, ".rds")
  )
  if (!file.exists(atac_file_sub)) {
    objects_atac <- readRDS(
      "../../data/BrainData/processed/GSE192774/GSE192774_ATAC_processed.rds"
    )
    objects_atac_sub <- subset(
      objects_atac, Sample == human_sample | Sample == chimp_sample
    )
    saveRDS(objects_atac_sub, atac_file_sub)
  } else {
    objects_atac_sub <- readRDS(atac_file_sub)
  }

  DefaultAssay(objects_atac_sub) <- "peaks"
  counts <- GetAssayData(objects_atac_sub, assay = "peaks", layer = "counts")
  peak_ids <- rownames(objects_atac_sub)

  log_message("Peak IDs: {.val {length(peak_ids)}}")

  parts <- strsplit(peak_ids, "[:-]+")
  len_ok <- lengths(parts) >= 3
  if (!all(len_ok)) {
    log_message(
      "Dropping {.val {sum(!len_ok)}} malformed peak IDs",
      message_type = "warning"
    )
    parts <- parts[len_ok]
    peak_ids <- peak_ids[len_ok]
  }
  idx <- !duplicated(peak_ids)
  parts <- parts[idx]
  peak_ids <- peak_ids[idx]
  mat <- do.call(rbind, parts)
  df <- data.frame(
    chr = mat[, 1],
    start = as.integer(mat[, 2]),
    end = as.integer(mat[, 3]),
    stringsAsFactors = FALSE
  )
  df <- df[!duplicated(paste(df$chr, df$start, df$end)), ]
  df <- df[order(df$chr, df$start), ]

  write.table(
    df[, c("chr", "start", "end")],
    file.path(res_dir, "union_peaks_hg38.bed"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  df$peak_id <- paste0(df$chr, ":", df$start, "-", df$end)
  write.csv(
    df,
    file.path(res_dir, "union_peaks_hg38.csv"),
    row.names = FALSE
  )

  log_message(
    "Completed union peaks for {.val {pair}}",
    message_type = "success"
  )
}
