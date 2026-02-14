source("code/functions/prepare_env.R")

library(GenomicRanges)
library(GenomicFeatures)
library(org.Hs.eg.db)

if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
}
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

sample_pairs <- list(
  list(human = "h4", chimp = "c4"),
  list(human = "h3", chimp = "c2"),
  list(human = "h1", chimp = "c1")
)

max_dist <- 500000L

log_message("Preparing TSS data from TxDb...")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
g <- genes(txdb)
g <- g[seqnames(g) %in% paste0("chr", c(1:22, "X", "Y"))]
sym <- suppressMessages(
  select(org.Hs.eg.db,
    keys = names(g),
    columns = "SYMBOL",
    keytype = "ENTREZID"
  )
)
g$symbol <- sym$SYMBOL[match(names(g), sym$ENTREZID)]
g <- g[!is.na(g$symbol)]
tss_pos <- ifelse(strand(g) == "+", start(g), end(g))
md <- data.frame(
  symbol = as.character(g$symbol),
  seqnames = as.character(seqnames(g)),
  tss = tss_pos,
  stringsAsFactors = FALSE
)
md <- md[!duplicated(md$symbol), ]
tss_gr <- GRanges(
  seqnames = md$seqnames,
  ranges = IRanges(md$tss, md$tss),
  symbol = md$symbol
)

for (pair in sample_pairs) {
  human_sample <- pair$human
  chimp_sample <- pair$chimp

  log_message(
    "Peak2Gene (distance <= {.val {max_dist/1000}} kb) for {.val {pair}}"
  )

  res_dir <- check_dir(
    paste0(
      "results/species_networks/", human_sample, "_", chimp_sample, "/atac"
    )
  )
  union_csv <- file.path(res_dir, "union_peaks_hg38.csv")

  if (!file.exists(union_csv)) {
    log_message(
      "Union peaks file not found: {.file {union_csv}}, skipping...",
      message_type = "warning"
    )
    next
  }

  peaks_df <- read.csv(union_csv, stringsAsFactors = FALSE)
  if (!all(c("chr", "start", "end", "peak_id") %in% colnames(peaks_df))) {
    if ("peak_id" %in% colnames(peaks_df)) {
      parts <- strsplit(peaks_df$peak_id, "[:-]+")
      len_ok <- lengths(parts) >= 3
      mat <- do.call(rbind, parts[len_ok])
      peaks_df <- data.frame(
        peak_id = peaks_df$peak_id[len_ok],
        chr = mat[, 1], start = as.integer(mat[, 2]),
        end = as.integer(mat[, 3]),
        stringsAsFactors = FALSE
      )
    }
  }

  peaks_gr <- GRanges(
    seqnames = peaks_df$chr,
    ranges = IRanges(peaks_df$start, peaks_df$end),
    peak_id = peaks_df$peak_id
  )

  tss_wide <- promoters(tss_gr, upstream = max_dist, downstream = max_dist)
  hits <- findOverlaps(peaks_gr, tss_wide)

  q <- queryHits(hits)
  s <- subjectHits(hits)
  tss_pos <- start(tss_gr[s])
  p_start <- start(peaks_gr[q])
  p_end <- end(peaks_gr[q])
  dist <- ifelse(tss_pos >= p_start & tss_pos <= p_end, 0,
    pmin(abs(p_start - tss_pos), abs(p_end - tss_pos))
  )

  out <- data.frame(
    peak_id = mcols(peaks_gr)$peak_id[q],
    gene = mcols(tss_gr)$symbol[s],
    distance = as.integer(dist),
    stringsAsFactors = FALSE
  )
  out <- out[out$distance <= max_dist, ]

  write.csv(
    out, file.path(res_dir, "peak2gene_500kb.csv"),
    row.names = FALSE,
    quote = FALSE
  )
  log_message(
    "Peak2Gene: {.val {nrow(out)}} peak-gene pairs (<= 500 kb) -> {.file {res_dir}}",
    message_type = "success"
  )
}
