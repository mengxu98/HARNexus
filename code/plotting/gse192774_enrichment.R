source("code/functions/prepare_env.R")

sample_pairs <- list(
  list(human = "h4", chimp = "c4"),
  list(human = "h3", chimp = "c2"),
  list(human = "h1", chimp = "c1")
)

human_color <- "#3271AE"
chimp_color <- "#D11A2D"

color_mapping <- c(
  "Human" = human_color,
  "Chimpanzee" = chimp_color
)

for (pair in sample_pairs) {
  human_sample <- pair$human
  chimp_sample <- pair$chimp

  log_message(
    "Loading data for {.val {pair}}..."
  )
  res_dir <- paste0(
    "results/species_networks/", human_sample, "_", chimp_sample, "/"
  )
  fig_dir <- check_dir(
    paste0(
      "figures/species_networks/", human_sample, "_", chimp_sample
    )
  )

  enrichment_by_species <- list()
  for (species in c("human", "chimp")) {
    sample <- ifelse(species == "human", human_sample, chimp_sample)
    object <- readRDS(
      file.path(
        res_dir,
        paste0(species, "_", sample, "_object.rds")
      )
    )

    genes <- read.csv(
      paste0(res_dir, "network_genes_by_celltype_", species, ".csv"),
      stringsAsFactors = FALSE
    )
    celltype_levels <- sort(as.character(unique(genes$CellType)))
    genes$CellType <- factor(
      genes$CellType,
      levels = celltype_levels
    )
    target_genes <- genes[!genes$is_TF, ]
    object$CellType <- factor(
      object$CellType,
      levels = celltype_levels
    )

    file_path <- file.path(
      fig_dir, paste0("enrichment_heatmap_", species, ".rds")
    )
    if (!file.exists(file_path)) {
      heatmap <- FeatureHeatmap(
        object,
        layer = "data",
        group.by = "CellType",
        nlabel = 0,
        features = target_genes$gene,
        feature_split = target_genes$CellType,
        group_palcolor = color_celltypes[celltype_levels],
        cell_annotation_palcolor = color_celltypes[celltype_levels],
        feature_split_palcolor = color_celltypes[celltype_levels],
        species = "Homo_sapiens",
        heatmap_palette = "viridis",
        db = "GO_BP",
        pvalueCutoff = 0.01,
        padjustCutoff = 0.05,
        Ensembl_version = 113,
        anno_terms = TRUE,
        terms_fontsize = 10,
        use_raster = TRUE,
        height = 6,
        width = 4
      )
      saveRDS(heatmap, file_path)
    } else {
      heatmap <- readRDS(file_path)
    }
    pdf(
      file.path(
        fig_dir,
        paste0("enrichment_heatmap_", species, ".pdf")
      ),
      width = 13.5, height = 7
    )
    print(heatmap$plot)
    dev.off()

    enrichment_by_species[[species]] <- heatmap$enrichment
  }

  pair_id <- paste0(human_sample, "_", chimp_sample)
  enrich_human <- enrichment_by_species$human$enrichment
  enrich_chimp <- enrichment_by_species$chimp$enrichment

  pvalue_cutoff <- 0.01
  padjust_cutoff <- 0.05
  enrich_human <- enrich_human[enrich_human$pvalue <= pvalue_cutoff & enrich_human$p.adjust <= padjust_cutoff, ]
  enrich_chimp <- enrich_chimp[enrich_chimp$pvalue <= pvalue_cutoff & enrich_chimp$p.adjust <= padjust_cutoff, ]


  df_human <- data.frame(
    sample_pair = pair_id,
    Species = "Human",
    CellType = as.character(enrich_human$Groups),
    ID = enrich_human$ID,
    Description = enrich_human$Description,
    is_brain_related = grepl(
      brain_pattern, enrich_human$Description,
      ignore.case = TRUE
    ),
    stringsAsFactors = FALSE
  )
  df_chimp <- data.frame(
    sample_pair = pair_id,
    Species = "Chimpanzee",
    CellType = as.character(enrich_chimp$Groups),
    ID = enrich_chimp$ID,
    Description = enrich_chimp$Description,
    is_brain_related = grepl(
      brain_pattern, enrich_chimp$Description,
      ignore.case = TRUE
    ),
    stringsAsFactors = FALSE
  )
  pathway_full <- rbind(df_human, df_chimp)

  by_cols <- c("sample_pair", "Species", "CellType")
  n_all <- aggregate(
    data.frame(n_pathways_all = pathway_full$ID),
    by = pathway_full[, by_cols],
    FUN = length
  )
  n_brain <- aggregate(
    data.frame(n_pathways_brain = pathway_full$is_brain_related),
    by = pathway_full[, by_cols],
    FUN = sum
  )
  pathway_compare <- merge(n_all, n_brain, by = by_cols)

  celltype_levels <- sort(as.character(unique(pathway_full$CellType)))
  pathway_full$CellType <- factor(
    pathway_full$CellType,
    levels = celltype_levels
  )
  pathway_compare$CellType <- factor(
    pathway_compare$CellType,
    levels = celltype_levels
  )
  pathway_compare$Species <- factor(
    pathway_compare$Species,
    levels = c("Human", "Chimpanzee")
  )
  pos_dodge <- position_dodge2(width = 0.8, preserve = "single")
  p_combined <- ggplot(pathway_compare, aes(x = CellType, fill = Species)) +
    geom_col(
      aes(y = n_pathways_all),
      position = pos_dodge,
      width = 0.75,
      alpha = 0.35,
      colour = "black",
      linewidth = 0.3
    ) +
    geom_col(
      aes(y = n_pathways_brain),
      position = pos_dodge,
      width = 0.75,
      alpha = 1,
      colour = "black",
      linewidth = 0.3
    ) +
    geom_point(
      data = data.frame(
        Pathway_type = c("Total", "Brain-related"),
        x = celltype_levels[1],
        y = 0
      ),
      aes(alpha = Pathway_type, x = x, y = y),
      size = 0,
      stroke = 0,
      key_glyph = draw_key_rect,
      inherit.aes = FALSE
    ) +
    geom_text(
      aes(y = n_pathways_all, label = n_pathways_all, group = Species),
      position = pos_dodge,
      vjust = -0.6,
      size = 2.5
    ) +
    geom_text(
      aes(y = n_pathways_brain, label = n_pathways_brain, group = Species),
      position = pos_dodge,
      vjust = 0,
      size = 2.5
    ) +
    scale_fill_manual(values = color_mapping, name = "Species") +
    scale_alpha_manual(
      values = c("Total" = 0.35, "Brain-related" = 1),
      name = "Pathway type",
      guide = guide_legend(override.aes = list(fill = "gray50"))
    ) +
    xlab("Celltype") +
    ylab("Number of pathways") +
    ylim(0, max(pathway_compare$n_pathways_all) * 1.05) +
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
    file.path(fig_dir, "enrichment_pathways_combined.pdf"),
    p_combined,
    width = 3.5,
    height = 3
  )

  log_message(
    "Completed heatmaps for {.val {pair}}",
    message_type = "success"
  )
}
