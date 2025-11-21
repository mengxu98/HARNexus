source("code/functions/utils.R")
source("code/functions/network_analysis.R")

results_dir <- check_dir("results/networks_pfc_ohnolog/")

ohnolog_data <- read.table(
  "data/genes_set/hsapiens.Pairs.Strict.2R_PMID31612943.txt",
  sep = "\t",
  header = TRUE
)

ohnolog_data <- ohnolog_data[!duplicated(ohnolog_data), ]
log_message(
  "Total ohnolog pairs: {.val {nrow(ohnolog_data)}}"
)

ohnolog_genes <- unique(
  c(ohnolog_data$Symbol1, ohnolog_data$Symbol2)
)
log_message(
  "Total unique ohnolog genes (symbols): {.val {length(ohnolog_genes)}}"
)

if (!file.exists(file.path(results_dir, "pfc_network.rds"))) {
  network_data <- read.csv(
    "data/networks/csv/network_data.csv"
  )
  pfc_network <- network_data[network_data$Region == "PFC", ]
  saveRDS(
    pfc_network, file.path(results_dir, "pfc_network.rds")
  )
} else {
  pfc_network <- readRDS(
    file.path(results_dir, "pfc_network.rds")
  )
}

seurat_objects <- list()
load("data/seurat_object_list/PFC_S8_seurat.Rdata")
seurat_objects[["PFC_S8"]] <- seurat
load("data/seurat_object_list/PFC_S9_seurat.Rdata")
seurat_objects[["PFC_S9"]] <- seurat


calculate_ohnolog_enrichment <- function(
    har_genes, ohnolog_genes, background_genes) {
  har_genes <- intersect(har_genes, background_genes)
  ohnolog_genes <- intersect(ohnolog_genes, background_genes)

  har_ohnolog <- length(
    intersect(har_genes, ohnolog_genes)
  )
  har_not_ohnolog <- length(
    setdiff(har_genes, ohnolog_genes)
  )
  non_har_ohnolog <- length(
    intersect(
      setdiff(background_genes, har_genes), ohnolog_genes
    )
  )
  non_har_not_ohnolog <- length(
    setdiff(background_genes, union(har_genes, ohnolog_genes))
  )

  contingency_table <- matrix(
    c(
      har_ohnolog, har_not_ohnolog,
      non_har_ohnolog, non_har_not_ohnolog
    ),
    nrow = 2, byrow = TRUE
  )

  fisher_result <- fisher.test(
    contingency_table,
    alternative = "greater"
  )

  har_prop <- har_ohnolog / length(har_genes)
  non_har_prop <- non_har_ohnolog / length(setdiff(background_genes, har_genes))
  fold_enrichment <- har_prop / non_har_prop

  list(
    har_ohnolog = har_ohnolog,
    har_total = length(har_genes),
    background_ohnolog = non_har_ohnolog,
    background_total = length(setdiff(background_genes, har_genes)),
    har_proportion = har_prop,
    background_proportion = non_har_prop,
    fold_enrichment = fold_enrichment,
    fisher_p = fisher_result$p.value,
    contingency_table = contingency_table
  )
}

permutation_test <- function(
    har_genes, ohnolog_genes, background_genes, n_perm = 1000) {
  observed_overlap <- sum(har_genes %in% ohnolog_genes)
  perm_overlaps <- numeric(n_perm)

  for (i in 1:n_perm) {
    random_genes <- sample(background_genes, length(har_genes))
    perm_overlaps[i] <- sum(random_genes %in% ohnolog_genes)
  }

  empirical_p <- (sum(perm_overlaps >= observed_overlap) + 1) / (n_perm + 1)

  list(
    observed_overlap = observed_overlap,
    perm_overlaps = perm_overlaps,
    empirical_p = empirical_p,
    mean_perm_overlap = mean(perm_overlaps),
    sd_perm_overlap = sd(perm_overlaps)
  )
}

gtf_data <- rtracklayer::import(
  "data/genes_set/Homo_sapiens.GRCh38.112.gtf.gz"
)
genes <- gtf_data[gtf_data$type == "gene"]
background_genes <- unique(genes$gene_name)
background_genes <- background_genes[!is.na(background_genes)]
log_message(
  "Background genes: {.val {length(background_genes)}}"
)


# stages_use <- c("S8", "S9")
stages_use <- sort(unique(pfc_network$Stage))
cell_types <- unique(pfc_network$CellType)

log_message(
  "Found {.val {length(stages_use)}} stages: {.val {stages_use}}"
)
log_message(
  "Found {.val {length(cell_types)}} cell types: {.val {cell_types}}"
)

cell_type_enrichment <- data.frame()
for (stage in stages_use) {
  stage_data <- pfc_network[pfc_network$Stage == stage, ]

  for (cell_type in cell_types) {
    cell_data <- stage_data[stage_data$CellType == cell_type, ]
    if (nrow(cell_data) == 0) next
    har_genes <- unique(cell_data$Target)

    enrich_result <- calculate_ohnolog_enrichment(
      har_genes, ohnolog_genes, background_genes
    )
    perm_result <- permutation_test(
      har_genes, ohnolog_genes, background_genes
    )

    cell_type_enrichment <- rbind(
      cell_type_enrichment, data.frame(
        Stage = stage,
        CellType = cell_type,
        HAR_genes_count = length(har_genes),
        fold_enrichment = enrich_result$fold_enrichment,
        p_value = enrich_result$fisher_p,
        empirical_p = perm_result$empirical_p,
        HAR_ohnolog_overlap = enrich_result$har_ohnolog,
        HAR_total = enrich_result$har_total,
        stringsAsFactors = FALSE
      )
    )
  }
}
print(cell_type_enrichment)
write.csv(
  cell_type_enrichment,
  file.path(results_dir, "cell_type_ohnolog_enrichment.csv"),
  row.names = FALSE
)

log_message("Starting pair-level analysis...")

classify_ohnolog_pairs <- function(ohnolog_data, har_genes) {
  pair <- data.frame(
    Symbol1 = ohnolog_data$Symbol1,
    Symbol2 = ohnolog_data$Symbol2,
    Gene1_HAR = ohnolog_data$Symbol1 %in% har_genes,
    Gene2_HAR = ohnolog_data$Symbol2 %in% har_genes,
    stringsAsFactors = FALSE
  )

  pair$classification <- ifelse(
    pair$Gene1_HAR & pair$Gene2_HAR, "Pair",
    ifelse(pair$Gene1_HAR | pair$Gene2_HAR, "Single", "None")
  )

  pair
}


pair_level_results <- data.frame()
for (stage in stages_use) {
  stage_data <- pfc_network[pfc_network$Stage == stage, ]

  for (cell_type in cell_types) {
    cell_data <- stage_data[stage_data$CellType == cell_type, ]
    if (nrow(cell_data) == 0) next

    har_genes <- unique(cell_data$Target)

    pair_classification <- classify_ohnolog_pairs(ohnolog_data, har_genes)
    pair_counts <- table(pair_classification$classification)

    pair_level_results <- rbind(
      pair_level_results,
      data.frame(
        Stage = stage,
        CellType = cell_type,
        HAR_genes_count = length(har_genes),
        Both_pairs = as.numeric(pair_counts["Pair"]),
        Single_pairs = as.numeric(pair_counts["Single"]),
        None_pairs = as.numeric(pair_counts["None"]),
        Total_pairs = sum(pair_counts),
        Pair_proportion = as.numeric(
          pair_counts["Pair"]
        ) / sum(pair_counts),
        Single_proportion = as.numeric(
          pair_counts["Single"]
        ) / sum(pair_counts),
        None_proportion = as.numeric(
          pair_counts["None"]
        ) / sum(pair_counts),
        stringsAsFactors = FALSE
      )
    )
  }
}
print(pair_level_results)
write.csv(
  pair_level_results,
  file.path(results_dir, "cell_type_ohnolog_pair_classification.csv"),
  row.names = FALSE
)

log_message("Analyzing ohnolog genes as network hub genes...")

library(igraph)
g <- graph_from_data_frame(
  pfc_network[, c("TF", "Target")],
  directed = TRUE
)

node_deg <- igraph::degree(g, mode = "all")
deg_names <- names(node_deg)

if (is.list(node_deg)) {
  node_deg <- node_deg[[1]]
}
node_deg <- as.numeric(node_deg)
names(node_deg) <- deg_names
degree_df <- data.frame(
  Gene = deg_names,
  Degree = node_deg,
  is_ohnolog = deg_names %in% ohnolog_genes,
  stringsAsFactors = FALSE
)

ohnolog_degrees <- degree_df$Degree[degree_df$is_ohnolog]
non_ohnolog_degrees <- degree_df$Degree[!degree_df$is_ohnolog]

degree_test <- wilcox.test(
  ohnolog_degrees, non_ohnolog_degrees,
  alternative = "greater"
)

high_degree_threshold <- as.numeric(
  quantile(degree_df$Degree, 0.9, names = FALSE, type = 7)
)
ohnolog_high_degree <- sum(ohnolog_degrees >= high_degree_threshold)
non_ohnolog_high_degree <- sum(non_ohnolog_degrees >= high_degree_threshold)

ohnolog_high_degree_prop <- if (length(ohnolog_degrees) > 0) {
  ohnolog_high_degree / length(ohnolog_degrees)
} else {
  NA_real_
}
non_ohnolog_high_degree_prop <- if (length(non_ohnolog_degrees) > 0) {
  non_ohnolog_high_degree / length(non_ohnolog_degrees)
} else {
  NA_real_
}
hub_fold_enrichment <- if (is.finite(non_ohnolog_high_degree_prop) && non_ohnolog_high_degree_prop > 0) {
  ohnolog_high_degree_prop / non_ohnolog_high_degree_prop
} else {
  NA_real_
}

log_message(
  sprintf(
    "High-degree genes (degree >= %d): Ohnolog %.1f%% vs Non-ohnolog %.1f%%, Fold enrichment: %s",
    high_degree_threshold,
    ohnolog_high_degree_prop * 100,
    non_ohnolog_high_degree_prop * 100,
    ifelse(is.na(hub_fold_enrichment), "NA", sprintf("%.2f", hub_fold_enrichment))
  )
)

hub_analysis_results <- data.frame(
  Analysis_Type = "Hub_Gene_Analysis",
  Metric = c(
    "Median_Degree_Ohnolog", "Median_Degree_NonOhnolog",
    "Wilcoxon_P_value", "High_Degree_Threshold",
    "Ohnolog_High_Degree_Prop", "NonOhnolog_High_Degree_Prop",
    "Hub_Fold_Enrichment"
  ),
  Value = c(
    median(ohnolog_degrees, na.rm = TRUE),
    median(non_ohnolog_degrees, na.rm = TRUE),
    degree_test$p.value,
    high_degree_threshold,
    ohnolog_high_degree_prop,
    non_ohnolog_high_degree_prop,
    hub_fold_enrichment
  ),
  Description = c(
    "Median degree of ohnolog genes",
    "Median degree of non-ohnolog genes",
    "Wilcoxon test p-value for degree comparison",
    "Threshold for high-degree genes (90th percentile)",
    "Proportion of ohnolog genes that are high-degree",
    "Proportion of non-ohnolog genes that are high-degree",
    "Fold enrichment of ohnolog genes among high-degree genes"
  ),
  stringsAsFactors = FALSE
)
write.csv(
  hub_analysis_results,
  file.path(results_dir, "ohnolog_hub_gene_analysis.csv"),
  row.names = FALSE
)


library(ggplot2)
library(patchwork)
library(dplyr)

plot_data <- data.frame(
  Degree = c(ohnolog_degrees, non_ohnolog_degrees),
  Group = c(
    rep("Ohnolog", length(ohnolog_degrees)),
    rep("Non-Ohnolog", length(non_ohnolog_degrees))
  )
)

p1 <- ggplot(
  plot_data, aes(x = Group, y = Degree, fill = Group)
) +
  geom_boxplot() +
  scale_y_log10() +
  scale_fill_manual(
    values = c("Ohnolog" = "#bf3232", "Non-Ohnolog" = "#117070")
  ) +
  theme_bw() +
  labs(
    x = "Gene Type",
    y = "Degree (log10 scale)",
    subtitle = sprintf("Wilcoxon test p = %.2e", degree_test$p.value)
  ) +
  theme(legend.position = "none")

enrichment_data <- data.frame(
  Group = c("Ohnolog", "Non-Ohnolog"),
  High_Degree_Prop = c(ohnolog_high_degree_prop, non_ohnolog_high_degree_prop),
  Count = c(ohnolog_high_degree, non_ohnolog_high_degree),
  Total = c(length(ohnolog_degrees), length(non_ohnolog_degrees))
)

p2 <- ggplot(
  enrichment_data,
  aes(x = Group, y = High_Degree_Prop, fill = Group)
) +
  geom_col() +
  geom_text(
    aes(
      label = sprintf("%d/%d\n(%.1f%%)", Count, Total, High_Degree_Prop * 100)
    ),
    vjust = -0.5, size = 3
  ) +
  scale_fill_manual(
    values = c("Ohnolog" = "#bf3232", "Non-Ohnolog" = "#117070")
  ) +
  ylim(0, max(enrichment_data$High_Degree_Prop, na.rm = TRUE) * 1.2) +
  theme_bw() +
  labs(
    x = "Gene Type",
    y = "Proportion of High-Degree Genes",
    subtitle = sprintf(
      "Fold enrichment: %s", sprintf("%.2f", hub_fold_enrichment)
    )
  ) +
  theme(legend.position = "none")

combined_plot <- p1 | p2

ggsave(
  file.path(results_dir, "fig.1-ohnolog_hub_gene_analysis.pdf"),
  combined_plot,
  width = 7, height = 4
)


library(pheatmap)
library(RColorBrewer)

# 1. Enrichment fold change heatmap
enrichment_matrix <- reshape2::dcast(
  cell_type_enrichment,
  Stage ~ CellType,
  value.var = "fold_enrichment"
)
rownames(enrichment_matrix) <- enrichment_matrix$Stage
enrichment_matrix <- enrichment_matrix[, -1]

enrichment_matrix <- as.matrix(enrichment_matrix)
enrichment_matrix[is.na(enrichment_matrix)] <- 1
enrichment_matrix[is.infinite(enrichment_matrix)] <- 1
enrichment_matrix[enrichment_matrix <= 0] <- 1

sig_matrix <- reshape2::dcast(
  cell_type_enrichment,
  Stage ~ CellType,
  value.var = "p_value"
)
rownames(sig_matrix) <- sig_matrix$Stage
sig_matrix <- sig_matrix[, -1]

sig_matrix <- as.matrix(sig_matrix)
sig_matrix[is.na(sig_matrix)] <- 1
sig_matrix <- sig_matrix < 0.05

colors <- colorRampPalette(c("white", "#bf3232"))(100)
pheatmap(
  enrichment_matrix,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  color = colors,
  annotation_row = NULL,
  annotation_col = NULL,
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color = "black",
  cellwidth = 25,
  cellheight = 25,
  main = "Ohnolog Enrichment in PFC Cell Types",
  filename = file.path(results_dir, "fig.2-cell_type_enrichment_heatmap.pdf"),
  width = 7, height = 7
)

# 2. Significant enrichment bar plot
sig_enrichment <- cell_type_enrichment[cell_type_enrichment$p_value < 0.05, ]
sig_enrichment <- sig_enrichment[order(sig_enrichment$fold_enrichment, decreasing = TRUE), ]

p_bar <- ggplot(
  sig_enrichment, aes(
    x = reorder(paste(Stage, CellType, sep = "_"), fold_enrichment),
    y = fold_enrichment, fill = CellType
  )
) +
  geom_col() +
  geom_text(aes(label = sprintf("%.2f", fold_enrichment)),
    hjust = 1, size = 2.2
  ) +
  coord_flip() +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  labs(
    title = "Ohnolog Enrichment",
    x = "Stage_CellType",
    y = "Fold Enrichment",
    fill = "Cell Type"
  ) +
  theme(
    title = element_text(size = 7),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    legend.position = "none"
  )

ggsave(
  file.path(results_dir, "fig.3-cell_type_enrichment_barplot.pdf"),
  p_bar,
  width = 3, height = 5
)

# 3. Enrichment dot plot
p_dot <- ggplot(
  cell_type_enrichment, aes(x = CellType, y = Stage)
) +
  geom_point(aes(size = HAR_ohnolog_overlap, color = -log10(p_value + 1e-300))) +
  scale_size_continuous(range = c(2, 8), name = "Overlap\nCount") +
  scale_color_gradient(
    low = "#0b617e", high = "#bf3232",
    name = "-log10(p-value)"
  ) +
  theme_bw() +
  labs(
    title = "Ohnolog Enrichment",
    x = "Cell Type",
    y = "Developmental Stage"
  ) +
  theme(
    title = element_text(size = 10),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8)
  )

ggsave(
  file.path(results_dir, "fig.4-cell_type_enrichment_dotplot.pdf"),
  p_dot,
  width = 5, height = 4
)




pair_prop_matrix <- reshape2::dcast(
  pair_level_results,
  Stage ~ CellType,
  value.var = "Pair_proportion"
)
rownames(pair_prop_matrix) <- pair_prop_matrix$Stage
pair_prop_matrix <- pair_prop_matrix[, -1]

single_prop_matrix <- reshape2::dcast(
  pair_level_results,
  Stage ~ CellType,
  value.var = "Single_proportion"
)
rownames(single_prop_matrix) <- single_prop_matrix$Stage
single_prop_matrix <- single_prop_matrix[, -1]

if (is.list(pair_prop_matrix)) {
  pair_prop_matrix <- as.data.frame(pair_prop_matrix)
}
if (is.list(single_prop_matrix)) {
  single_prop_matrix <- as.data.frame(single_prop_matrix)
}

pair_prop_matrix <- as.matrix(pair_prop_matrix)
single_prop_matrix <- as.matrix(single_prop_matrix)

pair_prop_matrix[is.na(pair_prop_matrix)] <- 0
single_prop_matrix[is.na(single_prop_matrix)] <- 0

library(ComplexHeatmap)
library(circlize)

max_pair <- max(pair_prop_matrix, na.rm = TRUE)
max_single <- max(single_prop_matrix, na.rm = TRUE)

ht1 <- Heatmap(
  as.matrix(pair_prop_matrix),
  name = "Both Pairs\nProportion",
  col = colorRamp2(c(0, max(max_pair, 0.01)), c("white", "#bf3232")),
  cluster_rows = FALSE,
  cluster_columns = ncol(pair_prop_matrix) > 1,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_title = "Developmental Stage",
  column_title = "Both Pairs Proportion",
  rect_gp = gpar(col = "grey30", lwd = 1)
)

ht2 <- Heatmap(
  as.matrix(single_prop_matrix),
  name = "Single Pairs\nProportion",
  col = colorRamp2(c(0, max(max_single, 0.01)), c("white", "#ffa500")),
  cluster_rows = FALSE,
  cluster_columns = ncol(single_prop_matrix) > 1,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_title = "Developmental Stage",
  column_title = "Single Pairs Proportion",
  rect_gp = gpar(col = "grey30", lwd = 1)
)

pdf(
  file.path(results_dir, "fig.5-ohnolog_pair_proportion_heatmap.pdf"),
  width = 6, height = 3.2
)
draw(ht1 + ht2, heatmap_legend_side = "right")
dev.off()



pair_summary <- pair_level_results %>%
  group_by(Stage) %>%
  summarise(
    Mean_Pair_Prop = mean(Pair_proportion, na.rm = TRUE),
    Mean_Single_Prop = mean(Single_proportion, na.rm = TRUE),
    Mean_None_Prop = mean(None_proportion, na.rm = TRUE),
    .groups = "drop"
  )

pair_summary_long <- reshape2::melt(
  pair_summary,
  id.vars = "Stage",
  variable.name = "Pair_Type",
  value.name = "Mean_Proportion"
)

pair_summary_long$Pair_Type <- factor(pair_summary_long$Pair_Type,
  levels = c("Mean_None_Prop", "Mean_Single_Prop", "Mean_Pair_Prop"),
  labels = c("None", "Single", "Both")
)

p_pair_summary <- ggplot(
  pair_summary_long, aes(x = Stage, y = Mean_Proportion, fill = Pair_Type)
) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c("None" = "#e0e0e0", "Single" = "#ffa500", "Both" = "#bf3232")) +
  theme_bw() +
  labs(
    x = "Developmental Stage",
    y = "Mean Proportion",
    fill = "Pair Type"
  ) +
  theme(legend.position = "right")

ggsave(
  file.path(results_dir, "fig.6-ohnolog_pair_summary.pdf"),
  p_pair_summary,
  width = 5, height = 3
)


all_cell_types <- sort(
  unique(c(pair_level_results$CellType, cell_type_enrichment$CellType))
)
all_stages <- sort(unique(pair_level_results$Stage))

complete_combinations <- expand.grid(
  Stage = all_stages, CellType = all_cell_types, stringsAsFactors = FALSE
)

pair_complete <- merge(
  complete_combinations, pair_level_results,
  by = c("Stage", "CellType"), all.x = TRUE
)
pair_complete[is.na(pair_complete)] <- 0

pair_long <- reshape2::melt(
  pair_complete[, c(
    "Stage", "CellType", "Pair_proportion",
    "Single_proportion", "None_proportion"
  )],
  id.vars = c("Stage", "CellType"),
  variable.name = "Pair_Type",
  value.name = "Proportion"
)

pair_long$Pair_Type <- factor(
  pair_long$Pair_Type,
  levels = c("None_proportion", "Single_proportion", "Pair_proportion"),
  labels = c("None", "Single", "Both")
)

pair_long$CellType <- factor(pair_long$CellType, levels = all_cell_types)

p_stacked <- ggplot(
  pair_long, aes(x = CellType, y = Proportion, fill = Pair_Type)
) +
  geom_col(position = "stack") +
  facet_wrap(~Stage, scales = "free_y") +
  scale_fill_manual(
    values = c("None" = "grey80", "Single" = "#ffa500", "Both" = "#bf3232")
  ) +
  scale_x_discrete(limits = all_cell_types) +
  theme_bw() +
  labs(
    title = "Ohnolog Pair Classification",
    x = "Cell Type",
    y = "Proportion",
    fill = "Pair Type"
  ) +
  theme(
    title = element_text(size = 10),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8)
  )

pair_count_long <- reshape2::melt(
  pair_complete[, c(
    "Stage", "CellType", "Both_pairs",
    "Single_pairs", "None_pairs"
  )],
  id.vars = c("Stage", "CellType"),
  variable.name = "Pair_Type",
  value.name = "Count"
)

pair_count_long$Pair_Type <- factor(
  pair_count_long$Pair_Type,
  levels = c("None_pairs", "Single_pairs", "Both_pairs"),
  labels = c("None", "Single", "Both")
)

pair_count_long$CellType <- factor(
  pair_count_long$CellType,
  levels = all_cell_types
)

p_count <- ggplot(
  pair_count_long, aes(x = CellType, y = Count, fill = Pair_Type)
) +
  geom_col(position = "dodge") +
  facet_wrap(~Stage, scales = "free_y") +
  scale_fill_manual(
    values = c("None" = "grey80", "Single" = "#ffa500", "Both" = "#bf3232")
  ) +
  scale_x_discrete(limits = all_cell_types) +
  theme_bw() +
  labs(
    title = "Ohnolog Pair Counts",
    x = "Cell Type",
    y = "Count",
    fill = "Pair Type"
  ) +
  theme(
    title = element_text(size = 10),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8)
  )

combined_pair_plot <- p_stacked + p_count +
  plot_layout(ncol = 2, guides = "collect")

ggsave(
  file.path(results_dir, "fig.7-ohnolog_pair_comparison.pdf"),
  combined_pair_plot,
  width = 12, height = 5
)
