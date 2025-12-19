# source("code/functions/prepare_env.R")

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

p_direction <- ggplot(
  grouped, aes(x = Direction, fill = Direction)
) +
  geom_bar(width = 0.5, color = "black") +
  scale_fill_manual(
    values = c(
      "gain_in_human" = "#e41a1c",
      "gain_in_chimp" = "#377eb8",
      "similar" = "gray60"
    ),
    labels = direction_labels
  ) +
  scale_x_discrete(labels = direction_labels) +
  theme_bw() +
  labs(
    x = "",
    y = "Number of HAR-TF Pairs",
    fill = "Direction"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(result_dir, "tf_gain_direction.pdf"),
  p_direction,
  width = 3.5,
  height = 3
)

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
p_col <- ggplot(
  top_gain[1:50, ],
  aes(x = TF_Name, y = N, fill = Direction)
) +
  geom_col(position = position_dodge()) +
  scale_fill_manual(
    values = c(
      "gain_in_human" = "#e41a1c",
      "gain_in_chimp" = "#377eb8"
    ),
    labels = direction_labels[c("gain_in_human", "gain_in_chimp")]
  ) +
  coord_flip() +
  theme_bw() +
  labs(
    x = "TFs",
    y = "Number of HARs with Motif Gain",
    fill = "Direction"
  )
ggsave(
  file.path(result_dir, "tf_gain_col.pdf"),
  p_col,
  width = 8,
  height = 6
)

tfs_human_gain <- top_gain[Direction == "gain_in_human", TF_Name]
tfs_chimp_gain <- top_gain[Direction == "gain_in_chimp", TF_Name]

# Create TF summary table with gain counts (Step 4)
tf_gain_summary <- grouped_gain[, .(
  gain_in_human = sum(Direction == "gain_in_human"),
  gain_in_chimp = sum(Direction == "gain_in_chimp")
), by = TF_Name]

# Calculate log2 ratio (Step 5)
tf_gain_summary[, log2_ratio := log2((gain_in_human + 1) / (gain_in_chimp + 1))]

# Save summary table
fwrite(
  tf_gain_summary,
  file.path(result_dir, "tf_gain_bias_summary.csv")
)

tf_gain_plot <- tf_gain_summary[gain_in_human + gain_in_chimp > 0]
tf_gain_plot <- setorder(tf_gain_plot, log2_ratio)

p_bias <- ggplot(
  tf_gain_plot,
  aes(
    x = forcats::fct_reorder(TF_Name, log2_ratio),
    y = log2_ratio, fill = log2_ratio > 0
  )
) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(
    values = c("FALSE" = "#377eb8", "TRUE" = "#e41a1c"),
    guide = "none"
  ) +
  coord_flip() +
  theme_bw() +
  labs(
    x = "TF",
    y = expression(log[2] * "(Human Gain / Chimp Gain)"),
    title = "TF Motif Evolution Bias in HARs"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  file.path(result_dir, "tf_gain_bias_barplot.pdf"),
  p_bias,
  width = 8,
  height = 20
)


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

p_venn <- ggVennDiagram(
  venn_list,
  category.names = names(venn_list),
  label_alpha = 0,
  label_color = "black",
  label_size = 4,
  set_color = colors_direction,
  edge_lty = "solid",
  edge_size = 1.5
) +
  scale_fill_gradient(
    low = "white",
    high = "white",
    guide = "none"
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    legend.position = "none"
  )

ggsave(
  file.path(result_dir, "tf_direction_venn.pdf"),
  p_venn,
  width = 4,
  height = 4
)

log_message("Creating heatmap...")

pivot_df <- dcast(
  grouped,
  HAR_ID ~ TF_Name,
  value.var = "Mean_Diff"
)

pivot_matrix <- as.matrix(pivot_df[, -1])
rownames(pivot_matrix) <- pivot_df$HAR_ID

tf_counts <- grouped[, .N, by = TF_Name]
valid_tfs <- tf_counts[N >= 0, TF_Name]
filtered_matrix <- pivot_matrix[, colnames(pivot_matrix) %in% valid_tfs]

filtered_matrix <- filtered_matrix[rowSums(!is.na(filtered_matrix)) > 0, ]

filtered_matrix[is.na(filtered_matrix)] <- 0

if (nrow(filtered_matrix) > 100) {
  variances <- apply(filtered_matrix, 1, var, na.rm = TRUE)
  top_100_hars <- names(sort(variances, decreasing = TRUE)[1:100])
  plot_matrix <- filtered_matrix[top_100_hars, ]
} else {
  plot_matrix <- filtered_matrix
}

plot_matrix[is.na(plot_matrix)] <- 0

har_variance <- apply(plot_matrix, 1, var, na.rm = TRUE)
har_mean_abs <- apply(plot_matrix, 1, function(x) mean(abs(x), na.rm = TRUE))
har_scores <- har_variance * har_mean_abs
top_20_hars <- names(
  sort(har_scores, decreasing = TRUE)[seq_len(min(20, length(har_scores)))]
)

tf_variance <- apply(plot_matrix, 2, var, na.rm = TRUE)
tf_mean_abs <- apply(plot_matrix, 2, function(x) mean(abs(x), na.rm = TRUE))
tf_scores <- tf_variance * tf_mean_abs
top_20_tfs <- names(
  sort(tf_scores, decreasing = TRUE)[seq_len(min(20, length(tf_scores)))]
)

col_fun <- colorRamp2(
  c(
    -max(abs(plot_matrix), na.rm = TRUE),
    0,
    max(abs(plot_matrix), na.rm = TRUE)
  ),
  c("#1966ad", "gray90", "#bb141a")
)

har_hist_data <- rowSums(abs(plot_matrix), na.rm = TRUE)
row_hist_anno <- rowAnnotation(
  Histogram = anno_barplot(
    har_hist_data,
    baseline = 0,
    gp = gpar(fill = "black", col = NA),
    width = unit(2, "cm"),
    axis_param = list(
      at = c(0, max(har_hist_data, na.rm = TRUE)),
      labels = c("0", sprintf("%.2f", max(har_hist_data, na.rm = TRUE)))
    )
  ),
  Top_20 = anno_mark(
    at = which(rownames(plot_matrix) %in% top_20_hars),
    labels = rownames(plot_matrix)[rownames(plot_matrix) %in% top_20_hars],
    side = "right",
    labels_gp = gpar(fontsize = 9, col = "black"),
    link_gp = gpar(lwd = 1, col = "#0c8b4e"),
    link_width = unit(5, "mm")
  ),
  width = unit(5, "cm")
)

tf_hist_data <- colSums(abs(plot_matrix), na.rm = TRUE)
col_hist_anno <- columnAnnotation(
  Histogram = anno_barplot(
    tf_hist_data,
    baseline = 0,
    gp = gpar(fill = "black", col = NA),
    height = unit(2, "cm"),
    axis_param = list(
      at = c(0, max(tf_hist_data, na.rm = TRUE)),
      labels = c("0", sprintf("%.2f", max(tf_hist_data, na.rm = TRUE)))
    )
  ),
  height = unit(2, "cm")
)
col_top_anno <- columnAnnotation(
  Top_20 = anno_mark(
    at = which(colnames(plot_matrix) %in% top_20_tfs),
    labels = colnames(plot_matrix)[colnames(plot_matrix) %in% top_20_tfs],
    side = "bottom",
    labels_gp = gpar(fontsize = 9, col = "black", fontface = "italic"),
    link_gp = gpar(lwd = 1, col = "#095c91"),
    link_height = unit(5, "mm")
  ),
  height = unit(5, "cm")
)

ht <- Heatmap(
  plot_matrix,
  name = "Mean_Diff",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  border = TRUE,
  # rect_gp = gpar(col = "gray", lwd = 0.05),
  column_title = "TF Binding Score Differences",
  row_title = "Top 100 HARs by Variance",
  heatmap_legend_param = list(
    title = "Mean\nDifference",
    at = c(
      -max(abs(plot_matrix), na.rm = TRUE),
      0,
      max(abs(plot_matrix), na.rm = TRUE)
    ),
    labels = c(
      sprintf("%.2f", -max(abs(plot_matrix), na.rm = TRUE)),
      "0",
      sprintf("%.2f", max(abs(plot_matrix), na.rm = TRUE))
    )
  ),
  na_col = "gray90",
  width = ncol(plot_matrix) * unit(1.1, "mm"),
  height = nrow(plot_matrix) * unit(1.1, "mm"),
  use_raster = FALSE,
  right_annotation = row_hist_anno,
  top_annotation = col_hist_anno,
  bottom_annotation = col_top_anno
)

pdf(
  file.path(result_dir, "motif_score_difference_heatmap_complex.pdf"),
  width = 15,
  height = 8
)
draw(ht)
dev.off()

log_message(
  sprintf(
    "Heatmap with top 20 annotations saved to %s",
    file.path(result_dir, "motif_score_difference_heatmap_complex.pdf")
  ),
  message_type = "success"
)

log_message("Motif processing complete!", message_type = "success")
