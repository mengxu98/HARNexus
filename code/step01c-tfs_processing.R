source("code/functions/prepare_env.R")

log_message("Processing HAR-TF pairs...")

colors_database <- c(
  "CIS_BP" = "#bf3553",
  "JASPAR" = "#FF6B35",
  "HOCOMOCO" = "#fecc11",
  "hTFTarget" = "#20a162",
  "Consensus" = "#2474b5"
)
for (species in c("human", "chimp")) {
  result_dir <- check_dir(
    file.path("results/step01-har_tf", species)
  )
  all_hits <- fread(file.path(result_dir, "combined_hits.csv"))
  har_info <- fread(file.path(result_dir, "har_coords.csv"))

  all_hits[, TF := tf_name]
  all_hits[, har := har_name]

  har_info[, Sequence := paste0(chrom, ":", start, "-", end)]
  all_hits <- merge(
    all_hits, har_info[, .(har, Sequence)],
    by.x = "har", by.y = "har", all.x = TRUE
  )

  cisbp_hits <- all_hits[database == "CISBP"]
  hocomoco_hits <- all_hits[database == "HOCOMOCO"]
  jaspar_hits <- all_hits[database == "JASPAR"]
  htftarget_hits <- all_hits[database == "HTFTARGET"]

  cisbp_summary <- cisbp_hits[, .(
    min_p = min(`p-value`, na.rm = TRUE),
    avg_score = mean(score, na.rm = TRUE),
    n_sites = .N,
    Sequence = first(Sequence),
    motif_name = first(motif_name),
    start = first(start),
    end = first(end),
    strand = first(strand)
  ), by = .(har, TF)]
  cisbp_summary[, source := "CIS_BP"]

  hocomoco_summary <- hocomoco_hits[, .(
    min_p = min(`p-value`, na.rm = TRUE),
    avg_score = mean(score, na.rm = TRUE),
    n_sites = .N,
    Sequence = first(Sequence),
    motif_name = first(motif_name),
    start = first(start),
    end = first(end),
    strand = first(strand)
  ), by = .(har, TF)]
  hocomoco_summary[, source := "HOCOMOCO"]

  jaspar_summary <- jaspar_hits[, .(
    min_p = min(`p-value`, na.rm = TRUE),
    avg_score = mean(score, na.rm = TRUE),
    n_sites = .N,
    Sequence = first(Sequence),
    motif_name = first(motif_name),
    start = first(start),
    end = first(end),
    strand = first(strand)
  ), by = .(har, TF)]
  jaspar_summary[, source := "JASPAR"]

  htftarget_summary <- htftarget_hits[, .(
    min_p = min(`p-value`, na.rm = TRUE),
    avg_score = mean(score, na.rm = TRUE),
    n_sites = .N,
    Sequence = first(Sequence),
    motif_name = first(motif_name),
    start = first(start),
    end = first(end),
    strand = first(strand)
  ), by = .(har, TF)]
  htftarget_summary[, source := "hTFTarget"]

  consensus_all <- merge(
    merge(
      merge(
        cisbp_summary[, .(
          har, TF,
          Sequence, start, end, strand,
          cisbp_score = avg_score, cisbp_p = min_p,
          cisbp_motif_name = motif_name
        )],
        jaspar_summary[, .(
          har, TF,
          jaspar_score = avg_score, jaspar_p = min_p,
          jaspar_motif_name = motif_name
        )],
        by = c("har", "TF"), all = FALSE
      ),
      hocomoco_summary[, .(
        har, TF,
        hocomoco_score = avg_score, hocomoco_p = min_p,
        hocomoco_motif_name = motif_name
      )],
      by = c("har", "TF"), all = FALSE
    ),
    htftarget_summary[, .(
      har, TF,
      htftarget_score = avg_score, htftarget_p = min_p,
      htftarget_motif_name = motif_name
    )],
    by = c("har", "TF"), all = FALSE
  )

  consensus_all[, jaspar_motif_name := gsub(" [^;]+", "", jaspar_motif_name)]

  consensus_all[, htftarget_motif_name := ifelse(
    grepl("Homo sapiens", htftarget_motif_name),
    htftarget_motif_name,
    NA_character_
  )]

  consensus_all[, htftarget_motif_name := gsub(
    "([A-Z0-9]+)%m-dataset-([0-9-]+)\\s+Homo\\s+sapiens",
    "\\1-dataset-\\2",
    htftarget_motif_name
  )]

  consensus_all[, htftarget_motif_name := sapply(
    strsplit(htftarget_motif_name, ";"),
    function(x) {
      formatted <- grep("-dataset-", x, value = TRUE)
      if (length(formatted) > 0) {
        return(paste(formatted, collapse = ";"))
      } else {
        simple_tf <- gsub("^([A-Z0-9]+)\\s+Homo\\s+sapiens.*$", "\\1", x)
        simple_tf <- simple_tf[simple_tf != x]
        if (length(simple_tf) > 0) {
          return(paste(unique(simple_tf), collapse = ";"))
        } else {
          return(NA_character_)
        }
      }
    }
  )]

  consensus_all[, htftarget_motif_name := gsub(
    "-dataset-",
    "-",
    htftarget_motif_name
  )]

  fwrite(
    consensus_all,
    file = file.path(result_dir, "har_tf_pairs_scores.csv"),
    sep = ","
  )
  fwrite(
    consensus_all[, .(
      har, Sequence, start, end, strand,
      cisbp_motif_name, jaspar_motif_name,
      hocomoco_motif_name, htftarget_motif_name, TF
    )],
    file = file.path(result_dir, "har_tf_pairs.csv"),
    sep = ","
  )


  cisbp_hars <- unique(cisbp_summary$har)
  jaspar_hars <- unique(jaspar_summary$har)
  hocomoco_hars <- unique(hocomoco_summary$har)
  htftarget_hars <- unique(htftarget_summary$har)
  all_hars <- unique(har_info$har)

  missing_hars_consensus <- setdiff(all_hars, unique(consensus_all$har))

  missing_hars_any <- setdiff(
    all_hars,
    union(union(union(cisbp_hars, jaspar_hars), hocomoco_hars), htftarget_hars)
  )

  missing_har_info_consensus <- har_info[har %in% missing_hars_consensus]

  fwrite(
    missing_har_info_consensus[, .(
      har, chrom, start, end, Sequence
    )],
    file = file.path(result_dir, "missing_hars_consensus.csv"),
    sep = ","
  )

  log_message("Total HARs in dataset: {.val {length(all_hars)}}")
  log_message("Consensus predictions: {.val {nrow(consensus_all)}}")
  log_message("Consensus HARs: {.val {uniqueN(consensus_all$har)}}")
  log_message("Consensus TFs: {.val {uniqueN(consensus_all$TF)}}")
  log_message("Missing consensus HARs: {.val {length(missing_hars_consensus)}}")

  cisbp_tfs <- unique(cisbp_summary$TF)
  jaspar_tfs <- unique(jaspar_summary$TF)
  hocomoco_tfs <- unique(hocomoco_summary$TF)
  htftarget_tfs <- unique(htftarget_summary$TF)
  consensus_tfs <- unique(consensus_all$TF)

  cisbp_pairs <- unique(
    paste(cisbp_summary$har, cisbp_summary$TF, sep = "::")
  )
  jaspar_pairs <- unique(
    paste(jaspar_summary$har, jaspar_summary$TF, sep = "::")
  )
  hocomoco_pairs <- unique(
    paste(hocomoco_summary$har, hocomoco_summary$TF, sep = "::")
  )
  htftarget_pairs <- unique(
    paste(htftarget_summary$har, htftarget_summary$TF, sep = "::")
  )
  consensus_pairs <- nrow(consensus_all)


  all_unique_tfs_consensus <- unique(
    c(cisbp_tfs, jaspar_tfs, hocomoco_tfs, htftarget_tfs)
  )
  upset_data_consensus <- data.frame(
    TF = all_unique_tfs_consensus,
    "CIS_BP" = as.integer(all_unique_tfs_consensus %in% cisbp_tfs),
    "JASPAR" = as.integer(all_unique_tfs_consensus %in% jaspar_tfs),
    "HOCOMOCO" = as.integer(all_unique_tfs_consensus %in% hocomoco_tfs),
    "hTFTarget" = as.integer(all_unique_tfs_consensus %in% htftarget_tfs),
    "Consensus" = as.integer(all_unique_tfs_consensus %in% consensus_tfs)
  )

  pdf(
    file.path(result_dir, "upset_tfs.pdf"),
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


  stats_summary <- data.table(
    Level = rep(c("TF-HAR Pairs", "HARs", "TFs"), each = 5),
    Method = rep(c("CIS_BP", "JASPAR", "HOCOMOCO", "hTFTarget", "Consensus"), 3),
    Count = c(
      length(cisbp_pairs),
      length(jaspar_pairs),
      length(hocomoco_pairs),
      length(htftarget_pairs),
      consensus_pairs,
      length(cisbp_hars),
      length(jaspar_hars),
      length(hocomoco_hars),
      length(htftarget_hars),
      uniqueN(consensus_all$har),
      length(cisbp_tfs),
      length(jaspar_tfs),
      length(hocomoco_tfs),
      length(htftarget_tfs),
      length(consensus_tfs)
    )
  )

  fwrite(
    stats_summary,
    file.path(result_dir, "statistics.csv")
  )

  multi_stats_vis <- data.table(
    Dimension = rep(c("TF-HAR Pairs", "HARs", "TFs"), each = 5),
    Method = rep(c("CIS_BP", "JASPAR", "HOCOMOCO", "hTFTarget", "Consensus"), 3),
    Count = c(
      length(cisbp_pairs), length(jaspar_pairs),
      length(hocomoco_pairs), length(htftarget_pairs), consensus_pairs,
      length(cisbp_hars), length(jaspar_hars),
      length(hocomoco_hars), length(htftarget_hars), uniqueN(consensus_all$har),
      length(cisbp_tfs), length(jaspar_tfs),
      length(hocomoco_tfs), length(htftarget_tfs), length(consensus_tfs)
    )
  )

  multi_stats_vis[, Dimension := factor(
    Dimension,
    levels = c("TF-HAR Pairs", "HARs", "TFs")
  )]
  multi_stats_vis[, Method := factor(
    Method,
    levels = c("CIS_BP", "JASPAR", "HOCOMOCO", "hTFTarget", "Consensus")
  )]

  plot_list <- lapply(
    c("TF-HAR Pairs", "HARs", "TFs"), function(i) {
      data <- multi_stats_vis[Dimension == i]
      ggplot(
        data, aes(x = Method, y = Count, fill = Method)
      ) +
        geom_col(width = 0.7) +
        geom_text(aes(label = Count), vjust = -0.5, size = 3) +
        scale_fill_manual(values = colors_database) +
        labs(title = paste0("Number of ", i), x = "", y = "Count") +
        theme_bw() +
        theme(
          legend.position = "bottom",
          axis.text.x = element_text(angle = 30, hjust = 1)
        ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
    }
  )

  combined_plot <- wrap_plots(plot_list) +
    plot_layout(ncol = 3, guides = "collect") &
    theme(legend.position = "bottom")

  ggsave(
    file.path(result_dir, "statistics.pdf"),
    combined_plot,
    width = 9, height = 3.5
  )

  venn_pairs_list <- list(
    "CIS_BP" = cisbp_pairs,
    "JASPAR" = jaspar_pairs,
    "HOCOMOCO" = hocomoco_pairs,
    "hTFTarget" = htftarget_pairs,
    "Consensus" = consensus_pairs
  )
  venn_pairs_plot <- ggVennDiagram(
    venn_pairs_list[1:4],
    category.names = c(
      names(venn_pairs_list[1:4])
    ),
    label_alpha = 0,
    label_color = "white",
    label_size = 5,
    set_color = colors_database[1:4],
    edge_lty = "solid",
    edge_size = 1.5
  ) +
    labs(title = "TF-HAR Pairs")

  ggsave(
    file.path(result_dir, "venn_tfhar_pairs_ggVennDiagram.pdf"),
    venn_pairs_plot,
    width = 7, height = 7
  )

  venn_pairs <- euler(
    venn_pairs_list
  )

  pdf(
    file.path(result_dir, "venn_tfhar_pairs.pdf"),
    width = 9, height = 9
  )
  print(
    plot(
      venn_pairs,
      quantities = list(type = c("counts", "percent"), cex = 1.0),
      fills = colors_database,
      edges = list(lwd = 1),
      labels = list(cex = 1.0, font = 2),
      main = "TF-HAR Pairs",
      cex.main = 1.2
    )
  )
  dev.off()

  venn_tfs_list <- list(
    "CIS_BP" = cisbp_tfs,
    "JASPAR" = jaspar_tfs,
    "HOCOMOCO" = hocomoco_tfs,
    "hTFTarget" = htftarget_tfs,
    "Consensus" = consensus_tfs
  )

  venn_tfs_plot <- ggVennDiagram(
    venn_tfs_list,
    category.names = c(
      names(venn_tfs_list)
    ),
    label_alpha = 0,
    label_color = "white",
    label_size = 5,
    set_color = colors_database,
    edge_lty = "solid",
    edge_size = 1.5
  ) +
    labs(title = "TFs")

  ggsave(
    file.path(result_dir, "venn_tfs_ggVennDiagram.pdf"),
    venn_tfs_plot,
    width = 8, height = 8
  )

  venn_tfs <- euler(venn_tfs_list)

  pdf(
    file.path(result_dir, "venn_tfs.pdf"),
    width = 8, height = 8
  )
  print(
    plot(
      venn_tfs,
      quantities = list(type = c("counts", "percent"), cex = 1.0),
      fills = colors_database,
      edges = list(lwd = 1),
      labels = list(cex = 1.0, font = 2),
      main = "Transcription Factors",
      cex.main = 1.2
    )
  )
  dev.off()


  top_n <- 10

  top_tfs <- head(consensus_all[, .N, by = TF][order(-N)], top_n)
  top_tfs[, TF := factor(TF, levels = rev(TF))]

  p_top_tfs <- ggplot(
    top_tfs, aes(x = TF, y = N)
  ) +
    geom_col(fill = "#0c8b4e", width = 0.5) +
    geom_text(aes(label = N), hjust = -0.2, size = 3) +
    coord_flip() +
    labs(
      title = paste0("Top ", top_n, " TFs"),
      x = "TFs",
      y = "Number of HARs"
    ) +
    theme_bw() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25)))

  top_hars <- head(consensus_all[, .N, by = har][order(-N)], top_n)
  top_hars[, har := factor(har, levels = rev(har))]

  p_top_hars <- ggplot(
    top_hars, aes(x = har, y = N)
  ) +
    geom_col(fill = "#095c91", width = 0.5) +
    geom_text(aes(label = N), hjust = -0.2, size = 3) +
    coord_flip() +
    labs(
      title = paste0("Top ", top_n, " HARs"),
      x = "HARs",
      y = "Number of TFs"
    ) +
    theme_bw() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25)))

  p_top <- p_top_hars + p_top_tfs +
    plot_layout(ncol = 2)

  ggsave(
    file.path(result_dir, paste0("top", top_n, "_hars_tfs.pdf")),
    p_top,
    width = 5, height = 2.5
  )
}


log_message("Creating combined statistics plot for human and chimp...")

human_stats <- fread(file.path("results/step01-har_tf/human/statistics.csv"))
human_stats[, Species := "Human"]

chimp_stats <- fread(file.path("results/step01-har_tf/chimp/statistics.csv"))
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
    data[Species == "Human", label_y := Count + max_count * 0.05]
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
        aes(x = x_offset, y = label_y, label = Count),
        hjust = 0.5, size = 2.7
      ) +
      geom_text(
        data = data[Species == "Chimp"],
        aes(x = x_offset, y = label_y, label = Count),
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
        axis.text.x = element_text(angle = 30, hjust = 1)
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.06)))

    species_dummy <- data.frame(
      Species = factor(c("Human", "Chimp"), levels = c("Human", "Chimp")),
      x = c(-Inf, -Inf),
      y = c(-Inf, -Inf)
    )

    p <- p +
      # pseudo points for species legend
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
  file.path(result_dir, "statistics_combined.pdf"),
  combined_all_plots,
  width = 10, height = 4
)
