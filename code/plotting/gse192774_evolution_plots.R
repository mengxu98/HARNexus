source("code/functions/prepare_env.R")

human_color <- "#3271AE"
chimp_color <- "#D11A2D"
color_mapping_species <- c(
  "Human" = human_color,
  "Chimpanzee" = chimp_color
)

peak_type_levels <- c(
  "Chromatin accessible",
  "Chromatin inaccessible"
)
color_mapping <- c(
  "Chromatin accessible" = "black",
  "Chromatin inaccessible" = "#CCCCCC"
)

color_mapping_traj <- c(
  "Regulatory rewiring" = "#006D87",
  "Regulatory innovation" = "#F9BD10",
  "Cis-regulatory activation" = "#0AA344",
  "Conserved" = "#5E7987"
)


sample_pairs <- list(
  list(human = "h4", chimp = "c4"),
  list(human = "h3", chimp = "c2"),
  list(human = "h1", chimp = "c1")
)

for (pair in sample_pairs) {
  log_message(
    "Creating comprehensive evolution plots for {.val {pair}}"
  )
  human_sample <- pair$human
  chimp_sample <- pair$chimp

  res_dir <- paste0(
    "results/species_networks/", human_sample, "_", chimp_sample, "/"
  )
  fig_dir <- check_dir(
    paste0("figures/species_networks/", human_sample, "_", chimp_sample)
  )

  evolution_data <- read.csv(
    file.path(
      res_dir, "evolution_daccre_for_target_genes.csv"
    ),
    stringsAsFactors = FALSE
  )
  evolution_data$Human_ATAC <- (evolution_data$n_Human_gained > 0) |
    (evolution_data$n_both_peaks > 0)
  evolution_data$Chimp_ATAC <- (evolution_data$n_Chimp_gained > 0) |
    (evolution_data$n_both_peaks > 0)

  in_human_net <- evolution_data$Target_type %in% c("Human-only", "Not biased")
  in_chimp_net <- evolution_data$Target_type %in% c("Chimp-only", "Not biased")
  n_total_human <- sum(in_human_net)
  n_accessible_human <- sum(evolution_data$Human_ATAC[in_human_net])
  n_total_chimp <- sum(in_chimp_net)
  n_accessible_chimp <- sum(evolution_data$Chimp_ATAC[in_chimp_net])
  atac_overall <- data.frame(
    Species = c("Human", "Chimpanzee"),
    n_total = c(n_total_human, n_total_chimp),
    n_accessible = c(n_accessible_human, n_accessible_chimp),
    stringsAsFactors = FALSE
  )
  atac_overall$Species <- factor(
    atac_overall$Species,
    levels = c("Human", "Chimpanzee")
  )
  pos_dodge <- position_dodge2(width = 0.8, preserve = "single")
  p_atac_by <- ggplot(atac_overall, aes(x = Species, fill = Species)) +
    geom_col(
      aes(y = n_total),
      position = pos_dodge,
      width = 0.75,
      alpha = 0.35,
      colour = "black",
      linewidth = 0.3
    ) +
    geom_col(
      aes(y = n_accessible),
      position = pos_dodge,
      width = 0.75,
      alpha = 1,
      colour = "black",
      linewidth = 0.3
    ) +
    geom_point(
      data = data.frame(
        Type = c("Total", "With chromatin accessibility"),
        x = atac_overall$Species[1],
        y = 0
      ),
      aes(alpha = Type, x = x, y = y),
      size = 0,
      stroke = 0,
      key_glyph = draw_key_rect,
      inherit.aes = FALSE
    ) +
    geom_text(
      aes(y = n_total * 0.98, label = n_total, group = Species),
      position = pos_dodge,
      vjust = -1,
      size = 2.5
    ) +
    geom_text(
      aes(y = n_accessible * 0.98, label = n_accessible, group = Species),
      position = pos_dodge,
      vjust = 1,
      size = 2.5
    ) +
    scale_fill_manual(values = color_mapping_species, name = "Species") +
    scale_alpha_manual(
      values = c("Total" = 0.35, "With chromatin accessibility" = 1),
      name = "Target genes",
      guide = guide_legend(override.aes = list(fill = "gray50"))
    ) +
    xlab("Species") +
    ylab("Number of target genes") +
    ylim(0, max(n_total_human, n_total_chimp) * 1.05) +
    theme_bw() +
    theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.key.size = unit(0.35, "cm"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  ggsave(
    file.path(fig_dir, "atac_type_distribution.pdf"),
    p_atac_by,
    width = 3,
    height = 2.5
  )

  celltype_levels <- sort(as.character(unique(evolution_data$CellType)))
  ct_has_chimp_data <- vapply(celltype_levels, function(ct) {
    any(evolution_data$Peak_type[evolution_data$CellType == ct] != "No_data")
  }, logical(1))
  atac_by_ct_list <- list()
  for (ct in celltype_levels) {
    ct_data <- evolution_data[evolution_data$CellType == ct, ]
    in_human_net_ct <- ct_data$Target_type %in% c("Human-only", "Not biased")
    in_chimp_net_ct <- ct_data$Target_type %in% c("Chimp-only", "Not biased")
    n_total_human_ct <- sum(in_human_net_ct)
    n_total_chimp_ct <- sum(in_chimp_net_ct)
    n_accessible_human_ct <- sum(ct_data$Human_ATAC[in_human_net_ct])
    n_accessible_chimp_ct <- sum(ct_data$Chimp_ATAC[in_chimp_net_ct])
    atac_by_ct_list[[length(atac_by_ct_list) + 1]] <- data.frame(
      CellType = ct,
      Species = "Human",
      n_total = n_total_human_ct,
      n_accessible = n_accessible_human_ct,
      stringsAsFactors = FALSE
    )
    if (ct_has_chimp_data[ct]) {
      atac_by_ct_list[[length(atac_by_ct_list) + 1]] <- data.frame(
        CellType = ct,
        Species = "Chimpanzee",
        n_total = n_total_chimp_ct,
        n_accessible = n_accessible_chimp_ct,
        stringsAsFactors = FALSE
      )
    }
  }
  atac_by_ct <- do.call(rbind, atac_by_ct_list)
  atac_by_ct$CellType <- factor(
    atac_by_ct$CellType,
    levels = celltype_levels
  )
  atac_by_ct$Species <- factor(
    atac_by_ct$Species,
    levels = c("Human", "Chimpanzee")
  )

  p_atac_by_ct <- ggplot(atac_by_ct, aes(x = CellType, fill = Species)) +
    geom_col(
      aes(y = n_total),
      position = pos_dodge,
      width = 0.75,
      alpha = 0.35,
      colour = "black",
      linewidth = 0.3
    ) +
    geom_col(
      aes(y = n_accessible),
      position = pos_dodge,
      width = 0.75,
      alpha = 1,
      colour = "black",
      linewidth = 0.3
    ) +
    geom_point(
      data = data.frame(
        Type = c("Total", "With chromatin accessibility"),
        x = celltype_levels[1],
        y = 0
      ),
      aes(alpha = Type, x = x, y = y),
      size = 0,
      stroke = 0,
      key_glyph = draw_key_rect,
      inherit.aes = FALSE
    ) +
    geom_text(
      aes(y = n_total, label = n_total, group = Species),
      position = pos_dodge,
      vjust = -0.6,
      size = 2.5
    ) +
    geom_text(
      aes(y = n_accessible * 0.98, label = n_accessible, group = Species),
      position = pos_dodge,
      vjust = 1,
      size = 2.5
    ) +
    scale_fill_manual(values = color_mapping_species, name = "Species") +
    scale_alpha_manual(
      values = c("Total" = 0.35, "With chromatin accessibility" = 1),
      name = "Target genes",
      guide = guide_legend(override.aes = list(fill = "gray50"))
    ) +
    xlab("Celltype") +
    ylab("Number of target genes") +
    ylim(0, max(atac_by_ct$n_total) * 1.05) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "right",
      legend.box = "vertical",
      legend.key.size = unit(0.35, "cm"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  ggsave(
    file.path(fig_dir, "atac_type_distribution_by_celltype.pdf"),
    p_atac_by_ct,
    width = 4.2,
    height = 3
  )

  trajectory_data <- read.csv(
    file.path(res_dir, "evolution_daccre_for_target_genes.csv"),
    stringsAsFactors = FALSE
  )
  table(trajectory_data$CellType) |> print()

  celltype_levels <- sort(as.character(unique(trajectory_data$CellType)))
  trajectory_data$CellType <- factor(
    trajectory_data$CellType,
    levels = celltype_levels
  )
  trajectory_data$Type <- factor(
    trajectory_data$Evolution_type,
    levels = names(color_mapping_traj)
  )

  trajectory_data_human <- trajectory_data[
    trajectory_data$Target_type %in% c("Human-only", "Not biased"),
  ]
  trajectory_data_human$CellType <- factor(
    as.character(trajectory_data_human$CellType),
    levels = sort(unique(as.character(trajectory_data_human$CellType)))
  )

  p_traj <- StatPlot(
    trajectory_data,
    stat.by = "Type",
    split.by = NULL,
    stat_type = "count",
    position = "dodge",
    label.size = 3,
    palcolor = color_mapping_traj,
    label = TRUE,
    bg.by = NULL
  )
  p_traj <- adjust_ggplot(p_traj)

  ggsave(
    file.path(fig_dir, "evolution_type_distribution.pdf"),
    p_traj,
    width = 4,
    height = 2.5
  )

  p_traj_by_ct <- StatPlot(
    trajectory_data_human,
    stat.by = "Type",
    plot_type = "bar",
    group.by = "CellType",
    # bg_alpha = 0.5,
    # bg_palcolor = color_celltypes[celltype_levels],
    # bg.by = "CellType",
    label.size = 3,
    palcolor = color_mapping_traj,
    stat_type = "count",
    position = "dodge",
    label = TRUE,
    xlab = "Celltype",
    ylab = "Number of target genes",
    theme_use = theme_bw
  )
  p_traj_by_ct <- adjust_ggplot(p_traj_by_ct)
  ggsave(
    file.path(fig_dir, "evolution_type_distribution_by_celltype.pdf"),
    p_traj_by_ct,
    width = 5,
    height = 3.5
  )
}
