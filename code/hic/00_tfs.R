tf_data <- read.csv("results/har_tf/human/har_tf_pairs_scores.csv")
tfs <- unique(tf_data$TF)
write.csv(
  tfs,
  file = "results/har_tf/tfs.csv",
  quote = FALSE,
  row.names = FALSE
)
