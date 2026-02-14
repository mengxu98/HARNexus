source("code/functions/prepare_env.R")

result_dir <- check_dir("results/har_tf")
result_dir_human <- check_dir(file.path(result_dir, "human"))
result_dir_chimp <- check_dir(file.path(result_dir, "chimp"))
genome_data_path <- "data/genome"

log_message("Starting HAR sequence extraction")

hars <- read.csv(
  file.path(genome_data_path, "HARs_PMID40011774.csv")
)
col_names <- c("chrom", "start", "end", "har")
har_human <- hars[, c(3, 4, 5, 2)]
colnames(har_human) <- col_names
har_chimp <- hars[, c(6, 7, 8, 2)]
colnames(har_chimp) <- col_names

write.csv(
  har_human,
  file.path(result_dir_human, "har_coords.csv"),
  row.names = FALSE
)

write.csv(
  har_chimp,
  file.path(result_dir_chimp, "har_coords.csv"),
  row.names = FALSE
)

gr_human <- GRanges(
  seqnames = har_human$chrom,
  ranges = IRanges(
    start = har_human$start, end = har_human$end
  )
)

seq_human <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr_human)
names(seq_human) <- har_human$har
writeXStringSet(
  seq_human,
  filepath = file.path(result_dir_human, "har_seqs_hg38.fa"),
  format = "fasta"
)

gr_chimp <- GRanges(
  seqnames = har_chimp$chrom,
  ranges = IRanges(
    start = har_chimp$start, end = har_chimp$end
  )
)
seq_chimp <- getSeq(BSgenome.Ptroglodytes.UCSC.panTro5, gr_chimp)
names(seq_chimp) <- har_chimp$har
writeXStringSet(
  seq_chimp,
  filepath = file.path(result_dir_chimp, "har_seqs_panTro5.fa"),
  format = "fasta"
)

log_message(
  "Extracted {.val {length(seq_human)}} and {.val {length(seq_chimp)}} human and chimp HAR sequences"
)
