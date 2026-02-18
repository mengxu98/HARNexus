source("code/functions/prepare_env.R")

colors_database <- c(
  "CIS_BP" = "#bf3553",
  "JASPAR" = "#FF6B35",
  "HOCOMOCO" = "#fecc11",
  "hTFTarget" = "#20a162",
  "Consensus" = "#2474b5"
)

har_color <- "#0C8B4E"
tf_color <- "#006D87"

for (species in c("human", "chimp")) {
  result_dir <- check_dir(file.path("results/har_tf", species))
  fig_dir <- check_dir(file.path("figures/har_tf", species))

  upset_data_consensus <- read.csv(
    file.path(result_dir, "upset_tfs_input.csv"),
    stringsAsFactors = FALSE
  )
  pdf(
    file.path(fig_dir, "upset_tfs.pdf"),
    width = 6, height = 5,
    onefile = FALSE
  )
  print(
    upset(
      upset_data_consensus,
      sets = c("CIS_BP", "JASPAR", "HOCOMOCO", "hTFTarget", "Consensus"),
      order.by = "freq",
      sets.bar.color = colors_database,
      main.bar.color = "gray30",
      matrix.color = "gray30",
      point.size = 5,
      line.size = 2,
      text.scale = c(1.5, 1.3, 1.3, 1.3, 1.3, 1.3),
      mainbar.y.label = "TF intersection count",
      sets.x.label = "TFs per database",
      mb.ratio = c(0.7, 0.3),
      queries = list(
        list(
          query = elements,
          params = list(
            "CIS_BP", "JASPAR",
            "HOCOMOCO", "hTFTarget", "Consensus"
          ),
          active = TRUE,
          color = "#378cc4",
          query.name = "Consensus"
        )
      )
    )
  )
  dev.off()


  stats_summary <- fread(
    file.path(result_dir, "statistics.csv")
  )
  stats_summary[, Dimension := factor(
    Level,
    levels = c("TF-HAR Pairs", "HARs", "TFs")
  )]
  stats_summary[, Method := factor(
    Method,
    levels = c("CIS_BP", "JASPAR", "HOCOMOCO", "hTFTarget", "Consensus")
  )]

  plot_list <- lapply(
    c("TF-HAR Pairs", "HARs", "TFs"), function(i) {
      data <- stats_summary[Dimension == i]
      ggplot(
        data, aes(x = Method, y = Count, fill = Method)
      ) +
        geom_rect(
          aes(
            xmin = as.numeric(Method) - 0.3,
            xmax = as.numeric(Method) + 0.3,
            ymin = 0,
            ymax = Count
          ),
          colour = "black",
          linewidth = 0.5
        ) +
        geom_text(aes(label = Count), vjust = -0.5, size = 3) +
        scale_fill_manual(values = colors_database) +
        labs(title = paste0("Number of ", i), x = "", y = "Count") +
        theme_bw() +
        theme(
          legend.position = "bottom",
          axis.text.x = element_text(angle = 30, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
    }
  )

  combined_plot <- wrap_plots(plot_list) +
    plot_layout(ncol = 3, guides = "collect") &
    theme(legend.position = "bottom")

  ggsave(
    file.path(fig_dir, "statistics.pdf"),
    combined_plot,
    width = 9, height = 3.5
  )

  venn_pairs_list <- readRDS(
    file.path(result_dir, "venn_tfhar_pairs_sets.rds")
  )
  venn_pairs_plot <- ggVennDiagram(
    venn_pairs_list[1:4],
    category.names = names(venn_pairs_list[1:4]),
    label_alpha = 0,
    label_color = "white",
    label_size = 5,
    set_color = colors_database[1:4],
    edge_lty = "solid",
    edge_size = 1.5
  ) +
    labs(title = "TF-HAR Pairs")

  ggsave(
    file.path(fig_dir, "venn_tfhar_pairs_ggVennDiagram.pdf"),
    venn_pairs_plot,
    width = 7, height = 7
  )

  venn_tfs_list <- readRDS(
    file.path(result_dir, "venn_tfs_sets.rds")
  )
  venn_tfs_plot <- ggVennDiagram(
    venn_tfs_list,
    category.names = names(venn_tfs_list),
    label_alpha = 0,
    label_color = "white",
    label_size = 5,
    set_color = colors_database,
    edge_lty = "solid",
    edge_size = 1.5
  ) +
    labs(title = "TFs")

  ggsave(
    file.path(fig_dir, "venn_tfs_ggVennDiagram.pdf"),
    venn_tfs_plot,
    width = 8, height = 8
  )
}

human_stats <- fread(
  "results/har_tf/human/statistics.csv"
)
human_stats[, Species := "Human"]

chimp_stats <- fread(
  "results/har_tf/chimp/statistics.csv"
)
chimp_stats[, Species := "Chimp"]

combined_stats <- rbind(human_stats, chimp_stats)
combined_stats[, Dimension := factor(
  Level,
  levels = c("TF-HAR Pairs", "HARs", "TFs")
)]
combined_stats[, Method := factor(
  Method,
  levels = c("CIS_BP", "JASPAR", "HOCOMOCO", "hTFTarget", "Consensus")
)]
combined_stats[, Species := factor(Species, levels = c("Human", "Chimp"))]

combined_plot_list <- lapply(
  c("TF-HAR Pairs", "HARs", "TFs"), function(dim) {
    data <- copy(combined_stats[Dimension == dim])
    method_levels <- levels(data$Method)
    data[, x_numeric := as.numeric(Method)]
    data[Species == "Human", x_offset := x_numeric - 0.2]
    data[Species == "Chimp", x_offset := x_numeric + 0.2]

    max_count <- max(data$Count, na.rm = TRUE)
    data[, label_y := Count + max_count * 0.05]
    data[Species == "Chimp", label_y := Count + max_count * 0.1]

    p <- ggplot(data, aes(x = x_offset, y = Count)) +
      geom_col(
        data = data[Species == "Human"],
        aes(fill = Method),
        width = 0.35
      ) +
      geom_col(
        data = data[Species == "Chimp"],
        aes(color = Method),
        fill = NA,
        width = 0.35
      ) +
      geom_text(
        data = data[Species == "Human"],
        aes(x = x_offset, y = label_y, label = format(Count, big.mark = ",")),
        hjust = 0.5, size = 2.7
      ) +
      geom_text(
        data = data[Species == "Chimp"],
        aes(x = x_offset, y = label_y, label = format(Count, big.mark = ",")),
        hjust = 0.5, size = 2.7
      ) +
      scale_fill_manual(values = colors_database, name = "Method") +
      scale_color_manual(values = colors_database, name = "Method") +
      scale_x_continuous(
        breaks = seq_along(method_levels),
        labels = method_levels
      ) +
      labs(
        title = paste0("Number of ", dim),
        x = "",
        y = "Count"
      ) +
      theme_bw() +
      theme(
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1)
      ) +
      scale_y_continuous(
        expand = expansion(mult = c(0, 0.06)),
        labels = function(x) format(x, scientific = FALSE, big.mark = ",")
      )

    species_dummy <- data.frame(
      Species = factor(c("Human", "Chimp"), levels = c("Human", "Chimp")),
      x = c(-Inf, -Inf),
      y = c(-Inf, -Inf)
    )

    p <- p +
      geom_point(
        data = species_dummy[species_dummy$Species == "Human", ],
        aes(x = x, y = y, shape = Species),
        fill = "black",
        color = "black",
        size = 0,
        alpha = 0,
        inherit.aes = FALSE
      ) +
      geom_point(
        data = species_dummy[species_dummy$Species == "Chimp", ],
        aes(x = x, y = y, shape = Species),
        fill = NA,
        color = "black",
        size = 0,
        alpha = 0,
        inherit.aes = FALSE
      ) +
      scale_shape_manual(
        values = c("Human" = 22, "Chimp" = 22),
        name = "Species"
      ) +
      guides(
        fill = guide_legend(title = "Method", order = 1),
        color = guide_legend(title = "Method", order = 1),
        shape = guide_legend(
          title = "Species",
          order = 2,
          override.aes = list(
            fill = c("black", NA),
            color = c("black", "black"),
            size = 8,
            shape = 22,
            alpha = 1
          )
        )
      )

    return(p)
  }
)

combined_all_plots <- wrap_plots(combined_plot_list) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  "figures/har_tf/statistics_combined.pdf",
  combined_all_plots,
  width = 11, height = 4
)
