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
  Mean_Diff > 0, "Gain in human",
  fifelse(Mean_Diff < 0, "Gain in chimpanzee", "No difference")
)]

fwrite(
  grouped,
  file.path(result_dir, "motif_score_aggregated.csv")
)

human_tfs <- unique(grouped[Direction == "Gain in human", TF_Name])
chimp_tfs <- unique(grouped[Direction == "Gain in chimpanzee", TF_Name])

summary_direction <- grouped[, .N, by = Direction]
log_message("Direction summary:")
print(summary_direction)

grouped_human <- grouped[Direction == "Gain in human"]
grouped_chimp <- grouped[Direction == "Gain in chimpanzee"]
grouped_similar <- grouped[Direction == "No difference"]

top_tf_gain_human <- grouped_human[, .N, by = TF_Name][order(-N)]
log_message("Top 10 TFs with gain in human:")
print(top_tf_gain_human[1:10])

top_tf_gain_chimp <- grouped_chimp[, .N, by = TF_Name][order(-N)]
log_message("Top 10 TFs with gain in chimp:")
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

grouped_gain <- grouped[Direction %in% c("Gain in human", "Gain in chimpanzee")]

top_gain <- grouped_gain[, .N, by = .(TF_Name, Direction)]

top_gain$TF_Name <- forcats::fct_reorder(
  top_gain$TF_Name, top_gain$N
)
top_gain <- setorder(top_gain, -N)
fwrite(
  top_gain,
  file.path(result_dir, "tf_gain_counts_by_direction.csv")
)

tfs_human_gain <- top_gain[Direction == "Gain in human", TF_Name]
tfs_chimp_gain <- top_gain[Direction == "Gain in chimpanzee", TF_Name]

tfs_gain_human <- unique(grouped[Direction == "Gain in human", TF_Name])
tfs_gain_chimp <- unique(grouped[Direction == "Gain in chimpanzee", TF_Name])
tfs_similar <- unique(grouped[Direction == "No difference", TF_Name])

log_message(
  sprintf(
    "TF counts - Human gain: %d, Chimp gain: %d, Similar: %d",
    length(tfs_gain_human),
    length(tfs_gain_chimp),
    length(tfs_similar)
  )
)

venn_list <- list(
  "Gain in human" = tfs_gain_human,
  "Gain in chimpanzee" = tfs_gain_chimp,
  "No difference" = tfs_similar
)

saveRDS(
  venn_list,
  file.path(result_dir, "tf_direction_sets.rds")
)

log_message("Motif processing complete!", message_type = "success")
