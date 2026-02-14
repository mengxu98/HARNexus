source("code/functions/prepare_env.R")

log_message("Plotting motif processing results...")


human_color <- "#3271AE"
chimp_color <- "#D11A2D"

har_color <- "#0C8B4E"
tf_color <- "#006D87"

result_dir <- check_dir("results/har_tf")
fig_dir <- check_dir("figures/har_tf")

grouped <- read.csv(
  file.path(result_dir, "motif_score_aggregated.csv"),
  stringsAsFactors = FALSE
)

direction_labels <- c(
  "gain_in_human" = "Gain in human",
  "gain_in_chimp" = "Gain in chimpanzee",
  "similar" = "No difference"
)

p_direction <- ggplot(
  grouped, aes(x = Direction, fill = Direction)
) +
  geom_bar(width = 0.5, color = "black") +
  geom_text(
    stat = "count",
    aes(label = after_stat(format(..count.., big.mark = ","))),
    vjust = -0.3,
    size = 3
  ) +
  scale_fill_manual(
    values = c(
      "gain_in_human" = human_color,
      "gain_in_chimp" = chimp_color,
      "similar" = "gray60"
    ),
    labels = direction_labels
  ) +
  scale_x_discrete(labels = direction_labels) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  ) +
  labs(
    x = "",
    y = "Number of HAR-TF pairs",
    fill = "Direction"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.08)),
    labels = function(x) format(x, scientific = FALSE, big.mark = ",")
  )

ggsave(
  file.path(fig_dir, "tf_gain_direction.pdf"),
  p_direction,
  width = 4,
  height = 3
)

venn_list <- readRDS(file.path(result_dir, "tf_direction_venn_sets.rds"))
colors_direction <- c(
  "Gain in human" = human_color,
  "Gain in chimpanzee" = chimp_color,
  "No difference" = "gray60"
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

label_mapping <- c(
  "Gain in Human" = "Gain in human",
  "Gain in Chimp" = "Gain in chimpanzee",
  "No Difference" = "No difference"
)

labels_to_move <- c(
  "Gain in human", "Gain in chimpanzee"
)

gb_venn <- ggplot_build(p_venn)
set_label_layer_idx <- NULL
set_label_layer_data <- NULL
for (i in seq_along(gb_venn$data)) {
  d <- gb_venn$data[[i]]
  if ("label" %in% names(d) && any(d$label %in% names(venn_list))) {
    set_label_layer_idx <- i
    set_label_layer_data <- d
    break
  }
}

if (!is.null(set_label_layer_idx) && !is.null(set_label_layer_data)) {
  set_labels_df <- unique(
    set_label_layer_data[set_label_layer_data$label %in% names(venn_list), c("x", "y", "label")]
  )

  set_labels_df$display_label <- sapply(set_labels_df$label, function(l) {
    if (l %in% names(label_mapping)) {
      label_mapping[[l]]
    } else {
      l
    }
  })

  center_x <- mean(set_labels_df$x, na.rm = TRUE)
  center_y <- mean(set_labels_df$y, na.rm = TRUE)

  should_move <- set_labels_df$label %in% labels_to_move |
    set_labels_df$display_label %in% c("Gain in human", "Gain in chimpanzee")

  set_labels_df$x <- ifelse(
    should_move,
    set_labels_df$x + (center_x - set_labels_df$x) * 0.4,
    set_labels_df$x
  )
  set_labels_df$y <- ifelse(
    should_move,
    set_labels_df$y,
    set_labels_df$y
  )

  p_venn$layers[[set_label_layer_idx]] <- NULL
  p_venn <- p_venn + geom_text(
    data = set_labels_df,
    aes(x = x, y = y, label = display_label),
    inherit.aes = FALSE,
    color = "black",
    size = 4
  )
}

ggsave(
  file.path(fig_dir, "venn_tf_direction.pdf"),
  p_venn,
  width = 3,
  height = 3
)

log_message("Creating heatmap...")
pivot_df <- dcast(
  grouped,
  HAR_ID ~ TF_Name,
  value.var = "Mean_Diff"
)
pivot_matrix <- as.matrix(pivot_df[, -1])
rownames(pivot_matrix) <- pivot_df$HAR_ID

filtered_matrix <- pivot_matrix[rowSums(!is.na(pivot_matrix)) > 0, , drop = FALSE]
filtered_matrix[is.na(filtered_matrix)] <- 0

if (nrow(filtered_matrix) > 100) {
  variances <- apply(filtered_matrix, 1, var, na.rm = TRUE)
  top_100_hars <- names(sort(variances, decreasing = TRUE)[1:100])
  plot_matrix <- filtered_matrix[top_100_hars, , drop = FALSE]
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
  c(chimp_color, "gray95", human_color)
)

har_hist_data <- rowSums(abs(plot_matrix), na.rm = TRUE)
row_hist_anno <- rowAnnotation(
  " " = anno_barplot(
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
    link_gp = gpar(lwd = 1, col = har_color),
    link_width = unit(5, "mm")
  ),
  width = unit(5, "cm")
)

tf_hist_data <- colSums(abs(plot_matrix), na.rm = TRUE)
col_hist_anno <- columnAnnotation(
  " " = anno_barplot(
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
bottom_anno <- columnAnnotation(
  Top_20 = anno_mark(
    at = which(colnames(plot_matrix) %in% top_20_tfs),
    labels = colnames(plot_matrix)[colnames(plot_matrix) %in% top_20_tfs],
    side = "bottom",
    labels_gp = gpar(fontsize = 9, col = "black", fontface = "italic"),
    link_gp = gpar(lwd = 1, col = tf_color),
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
  column_title = "TF binding score differences",
  row_title = "Top 100 HARs by variance",
  heatmap_legend_param = list(
    title = "Mean\ndifference",
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
  bottom_annotation = bottom_anno
)

pdf(
  file.path(fig_dir, "heatmap_tf_score_difference.pdf"),
  width = 15,
  height = 8
)
draw(ht)
dev.off()

log_message("Heatmap plotting complete", message_type = "success")
