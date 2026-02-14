source("code/functions/prepare_env.R")
source("code/functions/utils_network.R")

sample_pairs <- list(
  list(human = "h4", chimp = "c4"),
  list(human = "h3", chimp = "c2"),
  list(human = "h1", chimp = "c1")
)

tfs <- read.csv("results/har_tf/tfs.csv", stringsAsFactors = FALSE)[, 1]

for (pair in sample_pairs) {
  human_sample <- pair$human
  chimp_sample <- pair$chimp

  log_message(
    "Creating scatter plots for {.val {c(human_sample, chimp_sample)}}..."
  )
  res_dir <- paste0(
    "results/species_networks/", human_sample, "_", chimp_sample, "/"
  )
  fig_dir <- check_dir(
    paste0(
      "figures/species_networks/", human_sample, "_", chimp_sample
    )
  )
  csv_dir <- file.path(res_dir, "csv/")

  file_human <- file.path(
    res_dir, paste0("human_", human_sample, "_object.rds")
  )
  file_chimp <- file.path(
    res_dir, paste0("chimp_", chimp_sample, "_object.rds")
  )

  object_human <- readRDS(file_human)
  object_chimp <- readRDS(file_chimp)

  assay <- "RNA"
  layer <- "data"
  gene_sets_list_human <- get_celltype_specific_genes(
    object_human,
    assay = assay,
    layer = layer
  )
  gene_sets_list_chimp <- get_celltype_specific_genes(
    object_chimp,
    assay = assay,
    layer = layer
  )
  gene_sets_list <- purrr::map2(
    gene_sets_list_human, gene_sets_list_chimp, function(x, y) {
      intersect(x, y)
    }
  )

  plot_data_list <- list()

  human_files <- list.files(
    csv_dir,
    pattern = paste0("^human_", human_sample, "_.*\\.csv$"),
    full.names = TRUE
  )

  for (net_file in human_files) {
    ct <- gsub(paste0("^human_", human_sample, "_"), "", basename(net_file))
    ct <- gsub("\\.csv$", "", ct)

    net_df <- read.csv(net_file, stringsAsFactors = FALSE)
    if (nrow(net_df) == 0) next

    targets <- unique(net_df$target)
    targets <- targets[!targets %in% tfs]
    n_targets <- length(targets)

    n_cells <- sum(object_human$CellType == ct)

    n_specific_genes <- length(gene_sets_list[[ct]])

    plot_data_list[[length(plot_data_list) + 1]] <- data.frame(
      Species = "Human",
      CellType = ct,
      nCells = n_cells,
      nTargets = n_targets,
      nCellTypeSpecificGenes = n_specific_genes,
      stringsAsFactors = FALSE
    )
  }
  chimp_files <- list.files(
    csv_dir,
    pattern = paste0("^chimp_", chimp_sample, "_.*\\.csv$"),
    full.names = TRUE
  )

  for (net_file in chimp_files) {
    ct <- gsub(paste0("^chimp_", chimp_sample, "_"), "", basename(net_file))
    ct <- gsub("\\.csv$", "", ct)

    net_df <- read.csv(net_file, stringsAsFactors = FALSE)
    if (nrow(net_df) == 0) next

    targets <- unique(net_df$target)
    targets <- targets[!targets %in% tfs]
    n_targets <- length(targets)

    n_cells <- sum(object_chimp$CellType == ct)

    n_specific_genes <- length(gene_sets_list[[ct]])

    plot_data_list[[length(plot_data_list) + 1]] <- data.frame(
      Species = "Chimpanzee",
      CellType = ct,
      nCells = n_cells,
      nTargets = n_targets,
      nCellTypeSpecificGenes = n_specific_genes,
      stringsAsFactors = FALSE
    )
  }

  if (length(plot_data_list) == 0) {
    log_message(
      "No network data found for plotting",
      message_type = "warning"
    )
    next
  }

  plot_data <- do.call(rbind, plot_data_list)

  plot_data <- plot_data[plot_data$nTargets > 0, ]

  if (nrow(plot_data) == 0) {
    log_message(
      "No network data with targets found for plotting",
      message_type = "warning"
    )
    next
  }

  plot_data$CellType <- as.character(plot_data$CellType)

  cor_test1 <- cor.test(
    plot_data$nCells,
    plot_data$nTargets,
    method = "pearson"
  )
  cor_r1 <- round(cor_test1$estimate, 3)
  cor_p1 <- cor_test1$p.value
  cor_p1_label <- ifelse(cor_p1 < 0.001, "< 0.001", sprintf("%.3f", cor_p1))

  p1 <- ggplot(
    plot_data,
    aes(x = nCells, y = nTargets, color = CellType, shape = Species)
  ) +
    geom_point(size = 4) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      color = "black",
      linetype = "dashed",
      alpha = 0.2,
      inherit.aes = FALSE,
      aes(x = nCells, y = nTargets)
    ) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = paste0(
        "r = ", cor_r1, "\nP-value = ", cor_p1_label
      ),
      hjust = 1.1,
      vjust = 1.2,
      size = 4
    ) +
    scale_color_manual(
      values = color_celltypes,
      name = "Celltype",
      guide = guide_legend(override.aes = list(shape = 19))
    ) +
    scale_shape_manual(
      values = c("Human" = 19, "Chimpanzee" = 17),
      name = "Species",
      labels = c("Human" = "Human", "Chimpanzee" = "Chimpanzee")
    ) +
    labs(
      x = "Number of cells",
      y = "Number of target genes"
    ) +
    theme_bw() +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 13),
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 12),
      aspect.ratio = 1
    )

  cor_test2 <- cor.test(
    plot_data$nCellTypeSpecificGenes,
    plot_data$nTargets,
    method = "pearson"
  )
  cor_r2 <- round(cor_test2$estimate, 3)
  cor_p2 <- cor_test2$p.value
  cor_p2_label <- ifelse(cor_p2 < 0.001, "< 0.001", sprintf("%.3f", cor_p2))

  p2 <- ggplot(
    plot_data,
    aes(
      x = nCellTypeSpecificGenes,
      y = nTargets,
      color = CellType,
      shape = Species
    )
  ) +
    geom_point(size = 4) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      color = "black",
      linetype = "dashed",
      alpha = 0.2,
      inherit.aes = FALSE,
      aes(x = nCellTypeSpecificGenes, y = nTargets)
    ) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = paste0("r = ", cor_r2, "\nP-value = ", cor_p2_label),
      hjust = 1.1,
      vjust = 1.2,
      size = 4
    ) +
    scale_color_manual(
      values = color_celltypes,
      name = "Celltype",
      guide = guide_legend(override.aes = list(shape = 19))
    ) +
    scale_shape_manual(
      values = c("Human" = 19, "Chimpanzee" = 17),
      name = "Species",
      labels = c("Human" = "Human", "Chimpanzee" = "Chimpanzee")
    ) +
    labs(
      x = "Number of genes",
      y = ""
    ) +
    theme_bw() +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 13),
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 12),
      aspect.ratio = 1
    )

  p3 <- p1 + p2 + plot_layout(guides = "collect")
  ggsave(
    file.path(fig_dir, "scatter.pdf"),
    p3,
    width = 8.5,
    height = 4
  )
  log_message(
    "Scatter plots completed",
    message_type = "success"
  )
}
