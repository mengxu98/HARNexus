tf_data <- read.csv("results/har_tf/motif_score_aggregated.csv")
tfs <- unique(tf_data$TF_Name)
write.csv(
  tfs,
  file = "results/har_tf/tfs.csv",
  quote = FALSE,
  row.names = FALSE
)
