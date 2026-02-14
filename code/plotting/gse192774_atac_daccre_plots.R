source("code/functions/prepare_env.R")

sample_pairs <- list(
  list(human = "h4", chimp = "c4"),
  list(human = "h3", chimp = "c2"),
  list(human = "h1", chimp = "c1")
)

human_color <- "#3271AE"
chimp_color <- "#D11A2D"

fig_dir <- check_dir("figures/species_networks/")

p_list <- list()
for (pair in sample_pairs) {
  human_sample <- pair$human
  chimp_sample <- pair$chimp

  res_dir <- paste0(
    "results/species_networks/", human_sample, "_", chimp_sample, "/atac/"
  )

  log_message(
    "Creating DAcCRE visualizations for {.val {pair}}"
  )

  daccre_files <- list.files(
    res_dir,
    pattern = "^daccre_.*\\.csv$", full.names = TRUE
  )
  if (length(daccre_files) == 0) {
    log_message(
      "No DAcCRE files found in {.file {res_dir}}",
      message_type = "warning"
    )
    next
  }

  all_daccre <- list()
  celltypes <- character()

  for (f in daccre_files) {
    ct <- gsub("^daccre_", "", basename(f))
    ct <- gsub("\\.csv$", "", ct)
    ct <- gsub("_", " ", ct)
    if (ct == "Oligodendrocyte") {
      next
    }
    d <- read.csv(f, stringsAsFactors = FALSE)

    d$CellType <- ct
    all_daccre[[ct]] <- d
    celltypes <- c(celltypes, ct)
  }

  combined_daccre <- do.call(rbind, all_daccre)
  combined_daccre$CellType <- factor(combined_daccre$CellType)
  combined_daccre$peak_status <- factor(
    combined_daccre$peak_status,
    levels = c("Human-biased", "Chimp-biased", "Not-biased")
  )

  log_message("Loaded DAcCRE data for {.val {length(celltypes)}} cell types")

  status_summary <- as.data.frame(table(
    CellType = combined_daccre$CellType,
    Peak_type = combined_daccre$peak_status
  ))
  colnames(status_summary)[3] <- "Count"

  status_summary <- status_summary %>%
    group_by(CellType) %>%
    mutate(Proportion = Count / sum(Count)) %>%
    ungroup()
  diff_counts <- status_summary %>%
    filter(Peak_type %in% c("Human-biased", "Chimp-biased")) %>%
    dplyr::select(CellType, Peak_type, Count)

  p_diff_counts <- ggplot(
    diff_counts,
    aes(x = CellType, y = Count, fill = Peak_type)
  ) +
    geom_col(
      position = "dodge",
      width = 0.8,
      color = "black",
      linewidth = 0.3
    ) +
    geom_text(
      aes(label = Count),
      position = position_dodge(width = 0.9),
      vjust = -0.2,
      size = 2.5,
      color = "black"
    ) +
    scale_fill_manual(
      values = c(
        "Human-biased" = human_color,
        "Chimp-biased" = chimp_color
      ),
      name = "Peak type"
    ) +
    theme_bw() +
    ylim(0, max(diff_counts$Count) * 1.05) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    ) +
    labs(
      title = paste0("Human: ", human_sample, ", Chimpanzee: ", chimp_sample),
      x = "Celltype",
      y = "Number of differential peaks"
    )

  p_list[[paste0(human_sample, "_", chimp_sample)]] <- p_diff_counts
}

p_combined <- wrap_plots(p_list[1:2], ncol = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")
ggsave(
  file.path(fig_dir, "daccre_differential_counts.pdf"),
  p_combined,
  width = 8,
  height = 4
)
