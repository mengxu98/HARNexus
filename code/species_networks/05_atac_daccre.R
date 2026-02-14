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

  log_message(
    "DAcCRE (Human vs Chimp) per CellType for {.val {pair}}"
  )

  objects_atac_sub <- readRDS(
    file.path(
      res_dir,
      paste0(
        "GSE192774_ATAC_processed_", human_sample, "_", chimp_sample, ".rds"
      )
    )
  )

  counts <- GetAssayData(
    objects_atac_sub,
    assay = "peaks",
    layer = "counts"
  )

  meta <- objects_atac_sub@meta.data
  meta$CellType <- as.character(meta$CellType)
  meta$Species <- as.character(meta$Species)

  ct_use <- names(which(table(meta$CellType) >= 20))
  if (length(ct_use) == 0) stop("No cell type with >= 20 cells.")
  log_message("Cell types: {.val {paste(ct_use, collapse = ', ')}}")

  for (ct in ct_use) {
    idx <- meta$CellType == ct
    if (sum(idx) < 20) next
    sub_meta <- meta[idx, ]
    sub_counts <- counts[, idx, drop = FALSE]

    sp <- sub_meta$Species
    hu <- sp == "Human"
    ch <- sp == "Chimpanzee"
    if (sum(hu) < 5 || sum(ch) < 5) {
      log_message("  {.val {ct}}: too few Human/Chimp cells, skip")
      next
    }
    Y <- matrix(0, nrow = nrow(sub_counts), ncol = 2)
    Y[, 1] <- rowSums(sub_counts[, hu, drop = FALSE])
    Y[, 2] <- rowSums(sub_counts[, ch, drop = FALSE])
    rownames(Y) <- rownames(sub_counts)
    colnames(Y) <- c("Human", "Chimp")

    pseudocount <- 1
    log2fc <- log2((Y[, "Human"] + pseudocount) / (Y[, "Chimp"] + pseudocount))
    total_depth <- sum(Y)
    min_total <- max(20, round(total_depth / 50000))
    lfc_cut <- 0.5

    ok <- (Y[, "Human"] + Y[, "Chimp"]) >= min_total
    status <- rep("Not-biased", nrow(Y))
    status[ok & log2fc >= lfc_cut] <- "Human-biased"
    status[ok & log2fc <= -lfc_cut] <- "Chimp-biased"

    pid <- rownames(Y)
    if (length(pid) == 0) next
    parts <- strsplit(pid, "[:-]+")
    len_ok <- lengths(parts) >= 3

    mat <- do.call(rbind, parts[len_ok])
    chr <- mat[, 1]
    start <- as.integer(mat[, 2])
    end <- as.integer(mat[, 3])
    peak_id <- paste0(chr, ":", start, "-", end)
    out <- data.frame(
      peak_id = peak_id,
      chr = chr,
      start = start,
      end = end,
      Human_counts = Y[len_ok, "Human"],
      Chimp_counts = Y[len_ok, "Chimp"],
      Human_vs_Chimp_log2FC = log2fc[len_ok],
      peak_status = status[len_ok],
      stringsAsFactors = FALSE
    )
    if (sum(!len_ok) > 0) {
      out <- rbind(
        out, data.frame(
          peak_id = pid[!len_ok],
          chr = NA_character_, start = NA_integer_, end = NA_integer_,
          Human_counts = Y[!len_ok, "Human"],
          Chimp_counts = Y[!len_ok, "Chimp"],
          Human_vs_Chimp_log2FC = log2fc[!len_ok],
          peak_status = status[!len_ok],
          stringsAsFactors = FALSE
        )
      )
    }

    write.csv(
      out,
      file.path(res_dir, paste0("daccre_", gsub(" ", "_", ct), ".csv")),
      row.names = FALSE
    )
  }
}
