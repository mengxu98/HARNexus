source("code/functions/utils.R")
source("code/functions/network_analysis.R")

results_dir <- check_dir("results/networks_pfc_aging/")

colors_type <- c(
  "Aging" = "#3a83cc"
)

if (!file.exists(file.path(results_dir, "pfc_network.rds"))) {
  log_message("Processing data...")
  network_data <- read.csv("data/networks/csv/network_data.csv")
  pfc_network <- network_data[network_data$Region == "PFC", ]
  saveRDS(pfc_network, file.path(results_dir, "pfc_network.rds"))
} else {
  log_message("Loading data...")
  pfc_network <- readRDS(file.path(results_dir, "pfc_network.rds"))
}

if (!file.exists(file.path(results_dir, "enrichment_analysis.rds"))) {
  gtf_data <- import("data/genes_set/Homo_sapiens.GRCh38.112.gtf.gz")
  genes <- gtf_data[gtf_data$type == "gene"]
  background_genes <- unique(genes$gene_name)

  aging_gene <- fread("data/genes_set/genage_human.csv")
  aging_genes <- aging_gene$symbol

  ad_gene <- fread("data/genes_set/Alzheimer_Disease_gene.csv")
  ad_genes <- ad_gene$x

  pd_gene <- fread("data/genes_set/Parkinson_Disease_gene.csv")
  pd_genes <- pd_gene$x

  enrichment_data <- list(
    background_genes = background_genes,
    aging_genes = aging_genes,
    ad_genes = ad_genes,
    pd_genes = pd_genes
  )
  saveRDS(
    enrichment_data,
    file.path(results_dir, "enrichment_analysis.rds")
  )
} else {
  enrichment_data <- readRDS(
    file.path(results_dir, "enrichment_analysis.rds")
  )
  background_genes <- enrichment_data$background_genes
  aging_genes <- enrichment_data$aging_genes
}

stages <- sort(unique(pfc_network$Stage))
n_stages <- length(stages)
log_message(
  sprintf(
    "Found %d stages: %s",
    n_stages,
    paste(stages, collapse = ", ")
  )
)

log_message("Finding genes common to S8 and S9 but not in other stages...")

stage_genes_list <- list()
for (stage in stages) {
  stage_data <- pfc_network[pfc_network$Stage == stage, ]

  log_message(
    sprintf(
      "Stage %s: %d unique targets", stage, length(unique(stage_data$Target))
    )
  )

  stage_genes <- unique(stage_data$Target)
  stage_genes_list[[stage]] <- stage_genes
}

max_length <- max(sapply(stage_genes_list, length))
stage_genes_df <- data.frame(
  matrix(NA, nrow = max_length, ncol = length(stages))
)
colnames(stage_genes_df) <- stages

for (stage in stages) {
  genes <- stage_genes_list[[stage]]
  stage_genes_df[1:length(genes), stage] <- genes
}

write.xlsx(
  stage_genes_df,
  file.path(results_dir, "supp.file-stage_specific_genes.xlsx"),
  rowNames = FALSE
)


s8_genes <- unique(pfc_network[pfc_network$Stage == "S8", ]$Target)
s9_genes <- unique(pfc_network[pfc_network$Stage == "S9", ]$Target)

s8_s9_common_genes <- unique(c(s8_genes, s9_genes))
other_stages_genes <- unique(
  Reduce(
    union, lapply(
      setdiff(stages, c("S8", "S9")),
      function(stage) {
        unique(pfc_network[pfc_network$Stage == stage, ]$Target)
      }
    )
  )
)
s8_s9_unique_genes <- unique(
  setdiff(
    s8_s9_common_genes,
    other_stages_genes
  )
)

log_message(
  sprintf(
    "Found %d genes common to S8 and S9, %d unique to S8 and S9",
    length(s8_s9_common_genes),
    length(s8_s9_unique_genes)
  )
)

write.csv(
  s8_s9_common_genes,
  file.path(results_dir, "s8_s9_common_targets.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  s8_s9_unique_genes,
  file.path(results_dir, "s8_s9_unique_targets.csv"),
  row.names = FALSE,
  quote = FALSE
)


tf_lists <- lapply(
  stages, function(stage) {
    stage_data <- pfc_network[pfc_network$Stage == stage, ]
    unique(stage_data$TF)
  }
)
names(tf_lists) <- stages

target_lists <- lapply(
  stages, function(stage) {
    stage_data <- pfc_network[pfc_network$Stage == stage, ]
    unique(stage_data$Target)
  }
)
names(target_lists) <- stages

intersect_tfs <- Reduce(intersect, tf_lists)
intersect_targets <- Reduce(intersect, target_lists)

setdiff_targets_list <- lapply(
  stages, function(stage) {
    stage_data <- pfc_network[pfc_network$Stage == stage, ]
    stage_data_other <- pfc_network[pfc_network$Stage != stage, ]
    setdiff(stage_data$Target, stage_data_other$Target)
  }
)
names(setdiff_targets_list) <- stages
setdiff_targets <- Reduce(union, setdiff_targets_list)
write.csv(
  setdiff_targets,
  file.path(results_dir, "pfc_setdiff_targets.csv"),
  row.names = FALSE,
  quote = FALSE
)

all_targets <- Reduce(union, target_lists)
write.csv(
  all_targets,
  file.path(results_dir, "pfc_all_targets.csv"),
  row.names = FALSE,
  quote = FALSE
)

log_message(
  sprintf(
    "Found %d intersect TFs, and %d intersect targets...",
    length(intersect_tfs),
    length(intersect_targets)
  )
)

pfc_network_intersect <- pfc_network[pfc_network$Target %in% intersect_targets, ]
tf_lists_intersect <- lapply(
  stages, function(stage) {
    stage_data <- pfc_network_intersect[pfc_network_intersect$Stage == stage, ]
    unique(stage_data$TF)
  }
)
names(tf_lists_intersect) <- stages

write.csv(
  intersect_tfs,
  file.path(results_dir, "pfc_intersect_tfs.csv"),
  row.names = FALSE,
  quote = FALSE
)
write.csv(
  intersect_targets,
  file.path(results_dir, "pfc_intersect_targets.csv"),
  row.names = FALSE,
  quote = FALSE
)

stage_pairs <- list()
for (stage in stages) {
  stage_data <- pfc_network[pfc_network$Stage == stage, ]
  pairs <- paste(stage_data$TF, stage_data$Target, sep = "-")
  stage_pairs[[stage]] <- unique(pairs)
  log_message(
    sprintf(
      "Stage %s: %d unique TF-Target pairs",
      stage,
      length(unique(pairs))
    )
  )
}

jaccard_matrix <- matrix(0, nrow = n_stages, ncol = n_stages)
rownames(jaccard_matrix) <- stages
colnames(jaccard_matrix) <- stages

for (i in 1:n_stages) {
  for (j in 1:n_stages) {
    jaccard_matrix[i, j] <- calculate_jaccard(
      stage_pairs[[stages[i]]],
      stage_pairs[[stages[j]]]
    )
  }
}

pdf(
  file.path(results_dir, "fig.5a-jaccard_heatmap.pdf"),
  width = 4.5,
  height = 4,
  onefile = FALSE
)
plot_heatmap(jaccard_matrix)
dev.off()

jaccard_matrix_intersect <- calculate_jaccard_intersect(
  stage_pairs_list = stage_pairs,
  stages = stages,
  intersect_tfs = intersect_tfs,
  intersect_targets = intersect_targets
)

pdf(
  file.path(results_dir, "jaccard_heatmap_intersect.pdf"),
  width = 4.5,
  height = 4,
  onefile = FALSE
)
plot_heatmap(jaccard_matrix_intersect)
dev.off()

log_message("UpSet plots...")
pdf(
  file.path(results_dir, "fig.5b-upset_targets.pdf"),
  width = 8.5,
  height = 5.5,
  onefile = FALSE
)
print(upset_plot(target_lists))
dev.off()

pdf(
  file.path(results_dir, "upset_tfs.pdf"),
  width = 7.5,
  height = 6,
  onefile = FALSE
)
print(upset_plot(tf_lists))
dev.off()

pdf(
  file.path(results_dir, "upset_tfs_intersect.pdf"),
  width = 7.5,
  height = 6,
  onefile = FALSE
)
print(upset_plot(tf_lists_intersect))
dev.off()


tf_wordcloud_results <- plot_wordcloud(
  network_data = pfc_network_intersect,
  type = "TF",
  genes = intersect_tfs,
  min_freq = 3,
  max_words = 50
)

ggsave(
  file.path(results_dir, "tf_wordcloud.pdf"),
  tf_wordcloud_results$plot,
  width = 6,
  height = 6
)

target_wordcloud_results <- plot_wordcloud(
  network_data = pfc_network_intersect,
  type = "Target",
  genes = intersect_targets,
  min_freq = 3,
  max_words = 100
)

ggsave(
  file.path(results_dir, "target_wordcloud.pdf"),
  target_wordcloud_results$plot,
  width = 6,
  height = 6
)

tf_importance_results <- plot_importance(
  network_data = pfc_network_intersect,
  type = "TF",
  genes = intersect_tfs,
  top_n = 20
)

ggsave(
  file.path(results_dir, "tf_importance.pdf"),
  tf_importance_results$plot,
  width = 4,
  height = 3.5
)

target_importance_results <- plot_importance(
  network_data = pfc_network_intersect,
  type = "Target",
  genes = intersect_targets,
  top_n = 20
)

ggsave(
  file.path(results_dir, "target_importance.pdf"),
  target_importance_results$plot,
  width = 4,
  height = 3.5
)


log_message("Enrichment analysis...")

aging_enrichment <- calculate_enrichment(
  background_genes,
  s8_s9_common_genes,
  aging_genes
)

plot_data_combined <- data.frame(
  Group = factor(
    c("Background", "Aging"),
    levels = c("Background", "Aging")
  ),
  Fold_Enrichment = c(
    1, # Background is reference
    aging_enrichment$fold_enrichment
  ),
  p_value = c(
    1,
    aging_enrichment$p_value
  )
)

plot_data_combined$significance <- ifelse(
  plot_data_combined$p_value < 0.0001, "****",
  ifelse(plot_data_combined$p_value < 0.001, "***",
    ifelse(plot_data_combined$p_value < 0.01, "**",
      ifelse(plot_data_combined$p_value < 0.05, "*", "ns")
    )
  )
)
write.csv(
  plot_data_combined,
  file.path(results_dir, "enrichment_analysis.csv"),
  row.names = FALSE
)

plot_enrichment_all <- ggplot(
  plot_data_combined,
  aes(x = Group, y = Fold_Enrichment, fill = Group)
) +
  geom_bar(
    stat = "identity",
    width = 0.6,
    color = "black",
    size = 0.3
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "gray20"
  ) +
  theme_bw() +
  labs(
    y = "Fold enrichment",
    x = ""
  ) +
  ylim(0, max(plot_data_combined$Fold_Enrichment) * 1.05) +
  scale_fill_manual(
    values = c(
      "Background" = "gray30",
      colors_type
    )
  ) +
  geom_signif(
    comparisons = list(
      c("Background", "Aging")
    ),
    annotations = plot_data_combined$significance[-1],
    y_position = c(
      plot_data_combined$Fold_Enrichment[-1]
    ),
    map_signif_level = FALSE,
    tip_length = 0.01,
    size = 0.5,
    textsize = 3,
    vjust = 0.5,
    color = "black"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 0.8),
    legend.position = "right"
  )

ggsave(
  file.path(results_dir, "fig.5c-enrichment.pdf"),
  plot_enrichment_all,
  width = 3,
  height = 2.5
)


cell_types <- unique(pfc_network$CellType)
log_message("Calculating minimum edges needed to cover all genes...")
calculate_min_edges <- function(network_data, stages, cell_types) {
  min_edges <- list()

  for (stage in stages) {
    stage_data <- network_data[network_data$Stage == stage, ]
    stage_min_edges <- list()

    for (cell_type in cell_types) {
      log_message(
        sprintf(
          "Calculating minimum edges for %s in %s", cell_type, stage
        )
      )
      cell_data <- stage_data[stage_data$CellType == cell_type, ]
      if (nrow(cell_data) == 0) next

      all_targets <- unique(cell_data$Target)

      edges_needed <- 0
      covered_targets <- character(0)

      while (length(covered_targets) < length(all_targets)) {
        edges_needed <- edges_needed + 1
        current_edges <- head(cell_data, edges_needed)
        covered_targets <- unique(c(current_edges$Target, covered_targets))
      }

      stage_min_edges[[cell_type]] <- edges_needed
    }

    min_edges[[stage]] <- stage_min_edges
  }

  return(min_edges)
}

min_edges <- calculate_min_edges(pfc_network, stages_use, cell_types)

for (stage in stages_use) {
  log_message(sprintf("Stage %s:", stage))
  for (cell_type in names(min_edges[[stage]])) {
    log_message(
      sprintf(
        "  %s: %d edges needed", cell_type, min_edges[[stage]][[cell_type]]
      )
    )
  }
}

max_edges_needed <- max(unlist(min_edges))
log_message(
  sprintf(
    "Maximum edges needed across all stages and cell types: %d",
    max_edges_needed
  )
)

step_size <- 1000
edge_numbers <- seq(1000, 5000, by = step_size)
log_message(
  sprintf(
    "Using edge numbers from %d to %d with step size %d",
    edge_numbers[1], edge_numbers[length(edge_numbers)], step_size
  )
)

log_message("Analyzing enrichment across cell types for different stages...")
stages_use <- c("S8", "S9")
cell_types <- unique(pfc_network$CellType)

cell_type_enrichment_results <- list()

for (stage in stages_use) {
  stage_data <- pfc_network[pfc_network$Stage == stage, ]
  cell_type_results <- list()

  for (cell_type in cell_types) {
    cell_data <- stage_data[stage_data$CellType == cell_type, ]
    if (nrow(cell_data) == 0) next

    cell_results <- data.frame(
      n_edges = edge_numbers,
      aging_fold_enrichment = numeric(length(edge_numbers)),
      aging_p_value = numeric(length(edge_numbers))
    )

    for (i in seq_along(edge_numbers)) {
      n_edges <- edge_numbers[i]
      top_edges <- head(cell_data, n_edges)
      top_targets <- unique(top_edges$Target)

      aging_enrich <- calculate_enrichment(
        background_genes,
        top_targets,
        aging_genes
      )

      cell_results$aging_fold_enrichment[i] <- aging_enrich$fold_enrichment
      cell_results$aging_p_value[i] <- aging_enrich$p_value
    }

    cell_type_results[[cell_type]] <- cell_results
  }

  cell_type_enrichment_results[[stage]] <- cell_type_results
}

enrichment_summary <- data.frame()
for (stage in stages_use) {
  for (cell_type in names(cell_type_enrichment_results[[stage]])) {
    cell_results <- cell_type_enrichment_results[[stage]][[cell_type]]
    for (i in 1:nrow(cell_results)) {
      enrichment_summary <- rbind(
        enrichment_summary, data.frame(
          Stage = stage,
          CellType = cell_type,
          N_edges = cell_results$n_edges[i],
          Aging_fold_enrichment = cell_results$aging_fold_enrichment[i],
          Aging_p_value = cell_results$aging_p_value[i]
        )
      )
    }
  }
}

write.csv(
  enrichment_summary,
  file.path(results_dir, "cell_type_enrichment_summary.csv"),
  row.names = FALSE
)

plot_list <- list()
for (stage in stages_use) {
  stage_plots <- list()
  for (cell_type in names(cell_type_enrichment_results[[stage]])) {
    cell_results <- cell_type_enrichment_results[[stage]][[cell_type]]

    plot_data <- data.frame(
      n_edges = rep(cell_results$n_edges, 2),
      fold_enrichment = c(
        cell_results$aging_fold_enrichment
      ),
      p_value = c(
        cell_results$aging_p_value
      ),
      gene_set = rep(
        c("Aging"),
        each = nrow(cell_results)
      )
    )

    plot_data$significance <- ifelse(
      plot_data$p_value < 0.0001, "****",
      ifelse(plot_data$p_value < 0.001, "***",
        ifelse(plot_data$p_value < 0.01, "**",
          ifelse(plot_data$p_value < 0.05, "*", "ns")
        )
      )
    )

    p <- ggplot(
      plot_data,
      aes(x = n_edges, y = fold_enrichment, color = gene_set)
    ) +
      geom_line(size = 1) +
      geom_hline(
        yintercept = 1,
        linetype = "dashed",
        color = "gray20"
      ) +
      geom_point(size = 2) +
      geom_text(
        aes(label = significance),
        vjust = -0.2,
        size = 3,
        show.legend = FALSE
      ) +
      scale_color_manual(
        values = ifelse(stage == "S8", colors_type[1], colors_type[2]),
        name = "Gene set"
      ) +
      theme_bw() +
      labs(
        title = paste("Stage", stage, "-", cell_type),
        x = "Number of edges",
        y = "Fold enrichment"
      ) +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1)
      ) +
      scale_x_continuous(
        breaks = edge_numbers,
        labels = edge_numbers
      ) +
      ylim(0, max(plot_data$fold_enrichment) * 1.1)

    stage_plots[[cell_type]] <- p
    plot_list[[paste0(stage, "_", cell_type)]] <- p
  }
}

final_plot <- wrap_plots(plot_list, ncol = 5, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  file.path(results_dir, "s8_s9_cell_types_enrichment.pdf"),
  final_plot,
  width = 10,
  height = 6.5
)

log_message("Analyzing enrichment across cell types in S8 and S9 stages...")
stages_use <- c("S8", "S9")

cell_type_enrichment_summary <- data.frame()

for (stage in stages_use) {
  stage_data <- pfc_network[pfc_network$Stage == stage, ]

  for (cell_type in cell_types) {
    cell_data <- stage_data[stage_data$CellType == cell_type, ]
    if (nrow(cell_data) == 0) next

    targets <- unique(cell_data$Target)

    aging_enrich <- calculate_enrichment(
      background_genes,
      targets,
      aging_genes
    )

    cell_type_enrichment_summary <- rbind(
      cell_type_enrichment_summary, data.frame(
        Stage = stage,
        CellType = cell_type,
        Aging_fold_enrichment = aging_enrich$fold_enrichment,
        Aging_p_value = aging_enrich$p_value
      )
    )
  }
}

write.csv(
  cell_type_enrichment_summary,
  file.path(results_dir, "S8&S9_cell_type_enrichment_summary.csv"),
  row.names = FALSE
)

plot_data_simple <- cell_type_enrichment_summary %>%
  pivot_longer(
    cols = c(Aging_fold_enrichment),
    names_to = "Gene_Set",
    values_to = "Fold_Enrichment"
  ) %>%
  mutate(
    Gene_Set = ifelse(Gene_Set == "Aging_fold_enrichment", "Aging"),
    p_value = ifelse(
      Gene_Set == "Aging",
      Aging_p_value
    )
  )

plot_data_simple$significance <- ifelse(
  plot_data_simple$p_value < 0.0001, "****",
  ifelse(plot_data_simple$p_value < 0.001, "***",
    ifelse(plot_data_simple$p_value < 0.01, "**",
      ifelse(plot_data_simple$p_value < 0.05, "*", "ns")
    )
  )
)

colors_type_s8_s9 <- c(
  "S8" = "#3a83cc",
  "S9" = "#0b57a3"
)

plot_s8_s9_simple <- ggplot(
  plot_data_simple,
  aes(
    x = CellType,
    y = Fold_Enrichment,
    fill = Stage,
    group = interaction(Stage, Gene_Set)
  )
) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.9),
    width = 0.6,
    color = "gray20",
    size = 0.3
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "gray20"
  ) +
  geom_text(
    aes(label = significance),
    position = position_dodge(width = 0.9),
    vjust = -0.5,
    size = 3
  ) +
  facet_wrap(~Stage, ncol = 2) +
  scale_fill_manual(
    values = colors_type_s8_s9,
    name = "Stage"
  ) +
  theme_bw() +
  labs(
    x = "Cell type",
    y = "Fold enrichment"
  ) +
  theme(
    legend.position = "none"
  ) +
  ylim(0, max(plot_data_simple$Fold_Enrichment) * 1.15)

ggsave(
  file.path(results_dir, "S8&S9_cell_type_enrichment_summary.pdf"),
  plot_s8_s9_simple,
  width = 4.5,
  height = 2.2
)

plot_enrichment_combined <- plot_enrichment_all + plot_s8_s9_simple +
  plot_layout(ncol = 2, widths = c(0.2, 0.8)) +
  plot_annotation(tag_levels = "a")

ggsave(
  file.path(results_dir, "fig.5c-enrichment_combined.pdf"),
  plot_enrichment_combined,
  width = 7,
  height = 2.5
)


log_message("Performing enrichment analysis for S8 and S9 stages...")

s8_genes <- unique(pfc_network[pfc_network$Stage == "S8", ]$Target)
length(s8_genes)
s9_genes <- unique(pfc_network[pfc_network$Stage == "S9", ]$Target)
length(s9_genes)

s8_enrichment <- calculate_enrichment(
  background_genes,
  s8_genes,
  aging_genes
)

s9_enrichment <- calculate_enrichment(
  background_genes,
  s9_genes,
  aging_genes
)

plot_data_s8_s9 <- data.frame(
  Group = factor(
    c("Background", "S8", "S9"),
    levels = c("Background", "S8", "S9")
  ),
  Fold_Enrichment = c(
    1, # Background is reference
    s8_enrichment$fold_enrichment,
    s9_enrichment$fold_enrichment
  ),
  p_value = c(
    1,
    s8_enrichment$p_value,
    s9_enrichment$p_value
  )
)

plot_data_s8_s9$significance <- ifelse(
  plot_data_s8_s9$p_value < 0.0001, "****",
  ifelse(plot_data_s8_s9$p_value < 0.001, "***",
    ifelse(plot_data_s8_s9$p_value < 0.01, "**",
      ifelse(plot_data_s8_s9$p_value < 0.05, "*", "ns")
    )
  )
)

colors_type_s8_s92 <- c(
  "Background" = "gray50",
  "S8" = "#3a83cc",
  "S9" = "#0b57a3"
)

plot_enrichment_s8_s9 <- ggplot(
  plot_data_s8_s9,
  aes(x = Group, y = Fold_Enrichment, fill = Group)
) +
  geom_bar(
    stat = "identity",
    width = 0.6,
    color = "black",
    size = 0.3
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "gray20"
  ) +
  theme_bw() +
  labs(
    y = "Fold enrichment",
    x = ""
  ) +
  ylim(0, max(plot_data_s8_s9$Fold_Enrichment) * 1.1) +
  scale_fill_manual(
    values = colors_type_s8_s92,
    breaks = c("Background", "S8", "S9"),
    labels = c("Background", "S8", "S9")
  ) +
  geom_signif(
    comparisons = list(
      c("Background", "S8"),
      c("Background", "S9")
    ),
    annotations = plot_data_s8_s9$significance[-1],
    y_position = c(
      plot_data_s8_s9$Fold_Enrichment[-1]
    ),
    map_signif_level = FALSE,
    tip_length = 0.01,
    size = 0.5,
    textsize = 3,
    vjust = 0.5,
    color = "black"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "right"
  )

plot_enrichment_all_combined <- plot_enrichment_s8_s9 +
  plot_s8_s9_simple +
  plot_layout(ncol = 2, widths = c(0.2, 0.8)) +
  plot_annotation(tag_levels = "a")

ggsave(
  file.path(results_dir, "fig.5c-enrichment_all_combined.pdf"),
  plot_enrichment_all_combined,
  width = 8.5,
  height = 2.5
)

log_message("GO enrichment analysis...")
if (!file.exists(file.path(results_dir, "go_enrichment.rds"))) {
  go_bp <- enrichGO(
    gene = s8_s9_common_genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05
  )
  saveRDS(
    go_bp,
    file.path(results_dir, "go_enrichment.rds")
  )
} else {
  go_bp <- readRDS(
    file.path(results_dir, "go_enrichment.rds")
  )
}

go_bp_simp <- clusterProfiler::simplify(
  go_bp,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

go_bp_tree <- pairwise_termsim(
  go_bp,
  showCategory = 100
)
colors_go <- c(
  "#e95151",
  "#E69F00",
  "#56B4E9",
  "#009E73",
  "#F0E442",
  "#257a14"
)
cluster_n <- 5
plot_go_bp_tree <- treeplot(
  go_bp_tree,
  showCategory = 20,
  cluster.params = list(
    n = cluster_n,
    label_words_n = 3,
    color = colors_go[1:cluster_n],
    label_format = 5,
    hclust_method = "ward.D"
  ),
  clusterPanel.params = list(
    legend_n = 3
  ),
  geneClusterPanel = "heatMap"
)
ggsave(
  file.path(results_dir, "fig.5d-go_enrichment_s8_s9.pdf"),
  plot_go_bp_tree,
  width = 8,
  height = 6
)
write.csv(
  go_bp_simp,
  file.path(results_dir, "har_genes_go_enrichment_s8_s9.csv"),
  row.names = FALSE
)

log_message("GO enrichment analysis for S3-S7 stages...")

s3_s7_common_genes <- unique(
  Reduce(
    intersect,
    lapply(
      c("S3", "S4", "S5", "S6", "S7"),
      function(stage) {
        unique(pfc_network[pfc_network$Stage == stage, ]$Target)
      }
    )
  )
)

if (!file.exists(file.path(results_dir, "go_enrichment_s3_s7.rds"))) {
  go_bp_s3_s7 <- enrichGO(
    gene = s3_s7_common_genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05
  )
  saveRDS(
    go_bp_s3_s7,
    file.path(results_dir, "go_enrichment_s3_s7.rds")
  )
} else {
  go_bp_s3_s7 <- readRDS(
    file.path(results_dir, "go_enrichment_s3_s7.rds")
  )
}

go_bp_s3_s7_simp <- clusterProfiler::simplify(
  go_bp_s3_s7,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

go_bp_s3_s7_df <- as.data.frame(go_bp_s3_s7_simp)
go_bp_s3_s7_df <- go_bp_s3_s7_df[order(go_bp_s3_s7_df$p.adjust), ]
top_terms_s3_s7 <- head(go_bp_s3_s7_df, 20)


go_bp_s3_s7_tree <- pairwise_termsim(
  go_bp_s3_s7,
  showCategory = 100
)

plot_go_bp_s3_s7_tree <- treeplot(
  go_bp_s3_s7_tree,
  showCategory = 20,
  cluster.params = list(
    n = cluster_n,
    label_words_n = 3,
    color = colors_go[1:cluster_n],
    label_format = 5,
    hclust_method = "ward.D"
  ),
  clusterPanel.params = list(
    legend_n = 3
  ),
  geneClusterPanel = "heatMap"
)
ggsave(
  file.path(results_dir, "fig.5d-go_enrichment_s3_s7.pdf"),
  plot_go_bp_s3_s7_tree,
  width = 8,
  height = 6
)


log_message("Performing dynamic GO enrichment analysis...")
go_results_list <- list()
for (stage in stages) {
  network_stage <- pfc_network[pfc_network$Stage == stage, ]
  if (nrow(network_stage) > 1000) {
    network_stage <- network_stage[, c("TF", "Target", "Weight")]
    colnames(network_stage) <- c("regulator", "target", "weight")
    network_stage <- inferCSN::network_format(
      network_stage,
      abs_weight = FALSE
    )
    network_stage <- network_stage[1:1000, ]
  }
  stage_genes <- unique(network_stage$Target)

  if (length(stage_genes) == 0) {
    log_message(paste0("No stage-specific genes found for ", stage))
    next
  }

  log_message(
    paste0(
      "Processing ", stage, " stage-specific genes, ",
      length(stage_genes), " genes..."
    )
  )

  go_stage <- enrichGO(
    gene = stage_genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.2
  )

  go_results_list[[stage]] <- go_stage
}

go_results_list_top5 <- list()
for (stage in stages) {
  top_go <- head(go_results_list[[stage]]@result, 10)
  top_go$Stage <- stage

  go_results_list_top5[[stage]] <- top_go
}
dynamic_go_df <- do.call(rbind, go_results_list_top5)

go_matrix_count <- reshape2::dcast(
  dynamic_go_df,
  Description ~ Stage,
  value.var = "Count"
)
rownames(go_matrix_count) <- go_matrix_count$Description
go_matrix_count <- go_matrix_count[, -1]

go_matrix_count[is.na(go_matrix_count)] <- 1

go_matrix_count <- as.matrix(go_matrix_count)
p2_count <- pheatmap(
  go_matrix_count,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main = "GO enrichment of stage genes",
  fontsize_row = 8,
  fontsize_col = 8,
  color = colorRampPalette(c("#ffffff", "#1a5a9a"))(100),
  na_col = "white"
)
pdf(
  file.path(results_dir, "dynamic_go_heatmap_stage_specific.pdf"),
  width = 5,
  height = 5
)
print(p2_count)
dev.off()
