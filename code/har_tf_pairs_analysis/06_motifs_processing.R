source("code/functions/prepare_env.R")

log_message("Processing motif data...")

result_dir <- check_dir("results/har_tf")
df <- fread(file.path(result_dir, "motif_score_comparison.csv"))

grouped <- df[, .(
  Mean_Human_Score = mean(Human_Score, na.rm = TRUE),
  Mean_Chimp_Score = mean(Chimp_Score, na.rm = TRUE),
  Num_Motifs = .N
), by = .(HAR_ID, TF_Name)]

grouped[, Mean_Diff := Mean_Human_Score - Mean_Chimp_Score]
grouped[, Direction := fifelse(
  Mean_Diff > 0, "gain_in_human",
  fifelse(Mean_Diff < 0, "gain_in_chimp", "similar")
)]

direction_labels <- c(
  "gain_in_human" = "Gain in human",
  "gain_in_chimp" = "Gain in chimp",
  "similar" = "No difference"
)

fwrite(
  grouped,
  file.path(result_dir, "motif_score_aggregated.csv")
)

log_message("Aggregated results saved to {.file motif_score_aggregated.csv}")
human_tfs <- unique(grouped[Direction == "gain_in_human", TF_Name])
chimp_tfs <- unique(grouped[Direction == "gain_in_chimp", TF_Name])

summary_direction <- grouped[, .N, by = Direction]
log_message("Direction summary:", message_type = "info")
print(summary_direction)

grouped_human <- grouped[Direction == "gain_in_human"]
grouped_chimp <- grouped[Direction == "gain_in_chimp"]
grouped_similar <- grouped[Direction == "similar"]

top_tf_gain_human <- grouped_human[, .N, by = TF_Name][order(-N)]
log_message("Top 10 TFs with gain in human:", message_type = "info")
print(top_tf_gain_human[1:10])

top_tf_gain_chimp <- grouped_chimp[, .N, by = TF_Name][order(-N)]
log_message("Top 10 TFs with gain in chimp:", message_type = "info")
print(top_tf_gain_chimp[1:10])

har_gain_human <- grouped_human[, uniqueN(HAR_ID)]
har_gain_chimp <- grouped_chimp[, uniqueN(HAR_ID)]
har_similar <- grouped_similar[, uniqueN(HAR_ID)]

log_message(
  sprintf(
    "Unique HARs - Human gain: %d, Chimp gain: %d, Similar: %d",
    har_gain_human, har_gain_chimp, har_similar
  )
)

grouped_gain <- grouped[Direction %in% c("gain_in_human", "gain_in_chimp")]

top_gain <- grouped_gain[, .N, by = .(TF_Name, Direction)]

top_gain$TF_Name <- forcats::fct_reorder(
  top_gain$TF_Name, top_gain$N
)
top_gain <- setorder(top_gain, -N)
fwrite(
  top_gain,
  file.path(result_dir, "tf_gain_counts_by_direction.csv")
)

tfs_human_gain <- top_gain[Direction == "gain_in_human", TF_Name]
tfs_chimp_gain <- top_gain[Direction == "gain_in_chimp", TF_Name]

tf_gain_summary <- grouped_gain[, .(
  gain_in_human = sum(Direction == "gain_in_human"),
  gain_in_chimp = sum(Direction == "gain_in_chimp")
), by = TF_Name]

tf_gain_summary[, log2_ratio := log2((gain_in_human + 1) / (gain_in_chimp + 1))]

fwrite(
  tf_gain_summary,
  file.path(result_dir, "tf_gain_bias_summary.csv")
)

tf_gain_plot <- tf_gain_summary[gain_in_human + gain_in_chimp > 0]
tf_gain_plot <- setorder(tf_gain_plot, log2_ratio)


tfs_gain_human <- unique(grouped[Direction == "gain_in_human", TF_Name])
tfs_gain_chimp <- unique(grouped[Direction == "gain_in_chimp", TF_Name])
tfs_similar <- unique(grouped[Direction == "similar", TF_Name])

log_message(
  sprintf(
    "TF counts - Human gain: %d, Chimp gain: %d, Similar: %d",
    length(tfs_gain_human),
    length(tfs_gain_chimp),
    length(tfs_similar)
  )
)

venn_list <- list(
  "Gain in Human" = tfs_gain_human,
  "Gain in Chimp" = tfs_gain_chimp,
  "No Difference" = tfs_similar
)

colors_direction <- c(
  "Gain in Human" = "#e41a1c",
  "Gain in Chimp" = "#377eb8",
  "No Difference" = "gray60"
)
saveRDS(
  venn_list,
  file.path(result_dir, "tf_direction_venn_sets.rds")
)

log_message("Motif processing complete!", message_type = "success")
