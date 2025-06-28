source("code/functions/packages.R")

find_full_coverage_threshold <- function(
    weight_table,
    method_name) {
  all_target_genes <- unique(weight_table[[2]])
  total_possible_genes <- length(all_target_genes)

  log_message(method_name, " - Total unique target genes: ", total_possible_genes)

  thresholds <- seq(500, nrow(weight_table), by = 500)

  coverage_data <- data.frame(
    threshold = numeric(),
    genes_found = numeric(),
    percentage = numeric()
  )

  found_all_genes <- FALSE
  min_threshold_all_genes <- NA

  for (threshold in thresholds) {
    current_subset <- weight_table[1:threshold, ]
    current_genes <- length(unique(current_subset[[2]]))
    percentage <- current_genes / total_possible_genes

    coverage_data <- rbind(
      coverage_data,
      data.frame(
        threshold = threshold,
        genes_found = current_genes,
        percentage = percentage
      )
    )

    if (!found_all_genes && current_genes >= total_possible_genes) {
      found_all_genes <- TRUE
      min_threshold_all_genes <- threshold
      log_message(sprintf(
        "\n%s - Found all genes at threshold: %d",
        method_name, threshold
      ))
    }

    log_message(
      sprintf(
        "%s - Threshold: %d, Genes found: %d/%d (%.2f%%)",
        method_name, threshold,
        current_genes, total_possible_genes,
        percentage * 100
      )
    )
  }

  return(
    list(
      method = method_name,
      coverage_data = coverage_data,
      total_possible = total_possible_genes,
      min_threshold_all_genes = min_threshold_all_genes
    )
  )
}

coverage_trend_plot <- function(
    network_table,
    color_palette_method) {
  log_message("\nAnalyzing PPCOR:")
  ppcor_coverage <- find_full_coverage_threshold(
    network_table[["PPCOR"]],
    "PPCOR"
  )

  log_message("\nAnalyzing GENIE3:")
  genie3_coverage <- find_full_coverage_threshold(
    network_table[["GENIE3"]],
    "GENIE3"
  )

  log_message("\nAnalyzing LEAP:")
  leap_coverage <- find_full_coverage_threshold(
    network_table[["LEAP"]],
    "LEAP"
  )

  plot_data <- rbind(
    cbind(Method = "PPCOR", ppcor_coverage$coverage_data),
    cbind(Method = "GENIE3", genie3_coverage$coverage_data),
    cbind(Method = "LEAP", leap_coverage$coverage_data)
  )

  final_summary <- data.frame(
    Method = c("PPCOR", "GENIE3", "LEAP"),
    Total_Edges = c(
      nrow(network_table[["PPCOR"]]),
      nrow(network_table[["GENIE3"]]),
      nrow(network_table[["LEAP"]])
    ),
    Total_Genes = c(
      ppcor_coverage$total_possible,
      genie3_coverage$total_possible,
      leap_coverage$total_possible
    ),
    Min_edges_all_genes = c(
      ppcor_coverage$min_threshold_all_genes,
      genie3_coverage$min_threshold_all_genes,
      leap_coverage$min_threshold_all_genes
    ),
    Max_coverage = c(
      max(ppcor_coverage$coverage_data$percentage),
      max(genie3_coverage$coverage_data$percentage),
      max(leap_coverage$coverage_data$percentage)
    )
  )

  final_summary$Max_coverage_percent <- sprintf(
    "%.2f%%", final_summary$Max_coverage * 100
  )
  final_summary$Edges_ratio <- sprintf(
    "%.2f%%",
    final_summary$Min_edges_all_genes / final_summary$Total_Edges * 100
  )

  coverage_trend_plot <- ggplot(
    plot_data,
    aes(x = threshold, y = percentage, color = Method)
  ) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 0.7, alpha = 0.6) +
    geom_vline(
      data = final_summary,
      aes(xintercept = Min_edges_all_genes, color = Method),
      linetype = "dashed",
      alpha = 0.7,
      show.legend = FALSE
    ) +
    geom_text(
      data = final_summary,
      aes(
        x = Min_edges_all_genes,
        y = 0.5,
        label = sprintf("%d edges", Min_edges_all_genes),
        color = Method
      ),
      angle = 30,
      hjust = -0.1,
      vjust = 0,
      size = 2.5,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = color_palette_method[names(color_palette_method) %in% plot_data$Method]
    ) +
    scale_y_continuous(
      labels = scales::percent,
      limits = c(0, 1),
      breaks = seq(0, 1, 0.1)
    ) +
    scale_x_continuous(
      breaks = seq(0, max(plot_data$threshold), by = 10000)
    ) +
    labs(
      x = "Number of edges",
      y = "Genes coverage percentage",
      color = "Method"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1)
    )

  return(coverage_trend_plot)
}

calculate_intersection_trend <- function(
    network_table,
    tfdb_result,
    cbc_astro_hic,
    max_edges = 9000,
    step_size = 500) {
  results_df <- data.frame(
    threshold = numeric(),
    method = character(),
    total_genes_num = numeric(),
    intersect_genes_num = numeric(),
    ratio = numeric()
  )
  intersect_genes_list <- list()
  tfs_list <- list()
  all_genes_list <- list()

  network_table_infercsn <- network_table[["HARNexus"]]
  colnames(network_table_infercsn) <- c("TF", "gene", "weight")

  network_table_infercsn <- merge(
    network_table_infercsn,
    tfdb_result,
    by = "TF"
  )
  network_table_infercsn <- network_table_infercsn[, c(4, 5, 1, 2, 3)]

  compar_hic_infercsn <- merge(
    network_table_infercsn,
    cbc_astro_hic,
    by = "HAR"
  )

  filtered_df <- compar_hic_infercsn[compar_hic_infercsn$gene == compar_hic_infercsn$gene_hic, ]

  infercsn_total <- length(unique(compar_hic_infercsn$gene))
  infercsn_intersect <- length(unique(filtered_df$gene))
  infercsn_ratio <- infercsn_intersect / infercsn_total
  infercsn_row <- data.frame(
    threshold = 0,
    method = "HARNexus",
    total_genes_num = infercsn_total,
    intersect_genes_num = infercsn_intersect,
    ratio = infercsn_ratio
  )
  print(infercsn_row)
  tfs_list[["HARNexus"]][[as.character(0)]] <- as.character(
    unique(compar_hic_infercsn$TF)
  )
  all_genes_list[["HARNexus"]][[as.character(0)]] <- as.character(
    unique(compar_hic_infercsn$gene)
  )
  intersect_genes_list[["HARNexus"]][[as.character(0)]] <- as.character(
    unique(filtered_df$gene)
  )

  results_df <- rbind(results_df, infercsn_row)

  network_table_3 <- list(
    "PPCOR" = network_table[["PPCOR"]],
    "GENIE3" = network_table[["GENIE3"]],
    "LEAP" = network_table[["LEAP"]]
  )

  thresholds <- seq(step_size, max_edges, by = step_size)
  for (threshold in thresholds) {
    log_message("Calculating metrics for threshold: ", threshold)

    for (method_name in names(network_table_3)) {
      log_message("Processing ", method_name, " network...")
      weight_table <- network_table_3[[method_name]]

      current_threshold <- min(threshold, nrow(weight_table))
      weight_table_subset <- weight_table[1:current_threshold, ]

      colnames(weight_table_subset) <- c("TF", "gene", "weight")

      network <- merge(weight_table_subset, tfdb_result, by = "TF")
      network <- network[, c(4, 5, 1, 2, 3)]

      compar_hic <- merge(network, cbc_astro_hic, by = "HAR")

      filtered_df <- compar_hic[compar_hic[["gene"]] == compar_hic$gene_hic, ]

      total_genes_num <- length(unique(compar_hic[["gene"]]))
      intersect_genes <- unique(filtered_df[["gene"]])
      intersect_genes_list[[method_name]][[as.character(threshold)]] <- as.character(
        intersect_genes
      )
      tfs_list[[method_name]][[as.character(threshold)]] <- as.character(
        unique(compar_hic[["TF"]])
      )
      all_genes_list[[method_name]][[as.character(threshold)]] <- as.character(
        unique(compar_hic[["gene"]])
      )
      intersect_genes_num <- length(intersect_genes)

      ratio <- ifelse(
        total_genes_num > 0,
        intersect_genes_num / total_genes_num, 0
      )

      results_df <- rbind(
        results_df,
        data.frame(
          threshold = threshold,
          method = method_name,
          total_genes_num = total_genes_num,
          intersect_genes_num = intersect_genes_num,
          ratio = ratio
        )
      )

      log_message(
        sprintf(
          "    %s - threshold: %d, total genes: %d, intersect genes: %d, ratio: %.4f",
          method_name, threshold, total_genes_num, intersect_genes_num, ratio
        )
      )
    }
  }

  return(
    list(
      results = results_df,
      tfs_list = tfs_list,
      all_genes_list = all_genes_list,
      intersect_genes_list = intersect_genes_list
    )
  )
}

intersection_ratio_trend_plot <- function(
    intersection_results,
    color_palette_method) {
  intersection_ratio_trend_plot <- ggplot(
    intersection_results$results[-1, ],
    aes(x = threshold, y = ratio, color = method)
  ) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 0.7, alpha = 0.6) +
    labs(
      x = "Number of Edges",
      y = "Intersection Ratio",
      color = "Method"
    ) +
    theme_bw() +
    scale_color_manual(
      values = color_palette_method[names(color_palette_method) %in% intersection_results$results$method]
    ) +
    scale_y_continuous(
      labels = scales::percent,
      limits = c(0, 1),
      breaks = seq(0, 1, 0.1)
    ) +
    scale_x_continuous(
      breaks = seq(
        0,
        max(intersection_results$results$threshold),
        by = 1000
      )
    ) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 30, hjust = 1)
    )

  infercsn_ratio <- intersection_results$results$ratio[
    intersection_results$results$method == "HARNexus"
  ]
  intersection_ratio_trend_plot <- intersection_ratio_trend_plot +
    geom_hline(
      yintercept = infercsn_ratio,
      linetype = "solid",
      color = color_palette_method["HARNexus"]
    ) +
    annotate(
      "text",
      x = max(intersection_results$results$threshold),
      y = infercsn_ratio,
      label = sprintf("HARNexus (%.3f)", infercsn_ratio),
      hjust = 2,
      vjust = -0.5,
      color = color_palette_method["HARNexus"]
    ) +
    ylim(0, 0.3)

  return(intersection_ratio_trend_plot)
}

create_comparison_barplot <- function(
    results,
    thresholds = c(500, 1000, 2000),
    color_palette_method) {
  thresholds <- c(0, thresholds)
  filtered_data <- results$results[results$results$threshold %in% thresholds, ]

  filtered_data$threshold_label <- ifelse(
    filtered_data$method == "HARNexus",
    "HARNexus",
    paste0(filtered_data$threshold, " edges")
  )

  filtered_data$method <- factor(
    filtered_data$method,
    levels = c("GENIE3", "LEAP", "PPCOR", "HARNexus")
  )

  levels_threshold_label <- c(
    "HARNexus",
    paste0(unique(thresholds[-1]), " edges")
  )
  filtered_data$threshold_label <- factor(
    filtered_data$threshold_label,
    levels = levels_threshold_label
  )

  intersection_ratio_trend <- ggplot() +
    geom_bar(
      data = filtered_data,
      aes(x = method, y = ratio, fill = method),
      stat = "identity",
      width = 0.6
    ) +
    scale_fill_manual(
      values = color_palette_method,
      name = "Method"
    ) +
    labs(
      title = "",
      x = "Method",
      y = "Intersection ratio",
      fill = "Method"
    ) +
    geom_text(
      data = filtered_data,
      aes(
        x = method,
        y = ratio,
        label = sprintf(
          "%d/%d\n(%.1f%%)",
          intersect_genes_num,
          total_genes_num,
          ratio * 100
        )
      ),
      hjust = 1,
      vjust = 0.5,
      size = 5
    ) +
    coord_flip() +
    theme_bw() +
    theme(
      legend.position = "bottom",
      panel.grid.major.y = element_blank()
    ) +
    scale_y_continuous(
      limits = c(0, max(filtered_data$ratio) * 1.05),
      labels = scales::percent
    )

  return(intersection_ratio_trend)
}

extract_genes_list <- function(
    results,
    type = c("intersect", "all", "tfs", "setdiff"),
    threshold = 1000) {
  genes_list <- list()

  if (type == "intersect") {
    genes_list[["GENIE3"]] <- as.character(
      unique(results$intersect_genes_list$GENIE3[[as.character(threshold)]])
    )
    genes_list[["LEAP"]] <- as.character(
      unique(results$intersect_genes_list$LEAP[[as.character(threshold)]])
    )
    genes_list[["PPCOR"]] <- as.character(
      unique(results$intersect_genes_list$PPCOR[[as.character(threshold)]])
    )
    genes_list[["HARNexus"]] <- as.character(
      unique(results$intersect_genes_list$HARNexus[[as.character(0)]])
    )
  } else if (type == "all") {
    genes_list[["GENIE3"]] <- as.character(
      unique(results$all_genes_list$GENIE3[[as.character(threshold)]])
    )
    genes_list[["LEAP"]] <- as.character(
      unique(results$all_genes_list$LEAP[[as.character(threshold)]])
    )
    genes_list[["PPCOR"]] <- as.character(
      unique(results$all_genes_list$PPCOR[[as.character(threshold)]])
    )
    genes_list[["HARNexus"]] <- as.character(
      unique(results$all_genes_list$HARNexus[[as.character(0)]])
    )
  } else if (type == "tfs") {
    genes_list[["GENIE3"]] <- as.character(
      unique(results$tfs_list$GENIE3[[as.character(threshold)]])
    )
    genes_list[["LEAP"]] <- as.character(
      unique(results$tfs_list$LEAP[[as.character(threshold)]])
    )
    genes_list[["PPCOR"]] <- as.character(
      unique(results$tfs_list$PPCOR[[as.character(threshold)]])
    )
    genes_list[["HARNexus"]] <- as.character(
      unique(results$tfs_list$HARNexus[[as.character(0)]])
    )
  } else if (type == "setdiff") {
    genes_genie3_all <- as.character(
      unique(results$all_genes_list$GENIE3[[as.character(threshold)]])
    )
    genes_genie3_intersect <- as.character(
      unique(results$intersect_genes_list$GENIE3[[as.character(threshold)]])
    )
    genes_list[["GENIE3"]] <- as.character(
      setdiff(genes_genie3_all, genes_genie3_intersect)
    )
    genes_leap_all <- as.character(
      unique(results$all_genes_list$LEAP[[as.character(threshold)]])
    )
    genes_leap_intersect <- as.character(
      unique(results$intersect_genes_list$LEAP[[as.character(threshold)]])
    )
    genes_list[["LEAP"]] <- as.character(
      setdiff(genes_leap_all, genes_leap_intersect)
    )
    genes_ppcor_all <- as.character(
      unique(results$all_genes_list$PPCOR[[as.character(threshold)]])
    )
    genes_ppcor_intersect <- as.character(
      unique(results$intersect_genes_list$PPCOR[[as.character(threshold)]])
    )
    genes_list[["PPCOR"]] <- as.character(
      setdiff(genes_ppcor_all, genes_ppcor_intersect)
    )
    genes_infercsn_all <- as.character(
      unique(results$all_genes_list$HARNexus[[as.character(0)]])
    )
    genes_infercsn_intersect <- as.character(
      unique(results$intersect_genes_list$HARNexus[[as.character(0)]])
    )
    genes_list[["HARNexus"]] <- as.character(
      setdiff(genes_infercsn_all, genes_infercsn_intersect)
    )
  }

  return(genes_list)
}

create_venn_diagrams <- function(
    results,
    thresholds = 500,
    color_palette_method,
    return_upset = FALSE) {
  venn_plots <- list()

  for (threshold in thresholds) {
    log_message(paste("Creating Venn diagram for threshold", threshold, "..."))

    method_genes <- extract_genes_list(
      results,
      type = "intersect",
      threshold = threshold
    )
    method_genes <- list(
      HARNexus = method_genes$HARNexus,
      GENIE3 = method_genes$GENIE3,
      PPCOR = method_genes$PPCOR,
      LEAP = method_genes$LEAP
    )
    color_palette_method <- c(
      "HARNexus" = color_palette_method["HARNexus"],
      "GENIE3" = color_palette_method["GENIE3"],
      "PPCOR" = color_palette_method["PPCOR"],
      "LEAP" = color_palette_method["LEAP"]
    )

    venn_plot <- ggVennDiagram::ggVennDiagram(
      method_genes,
      category.names = names(method_genes),
      label_alpha = 0,
      label_color = "white",
      label_size = 5,
      set_color = color_palette_method,
      edge_lty = "solid",
      edge_size = 1.5
    )

    if (return_upset) {
      venn_plot <- ggVennDiagram::plot_upset(
        ggVennDiagram::Venn(method_genes),
        nintersects = NULL,
        order.intersect.by = c("size", "name", "none"),
        order.set.by = c("size", "name", "none"),
        relative_height = 3,
        relative_width = 0.3,
        top.bar.color = "#238ace",
        top.bar.y.label = NULL,
        top.bar.show.numbers = TRUE,
        top.bar.numbers.size = 5,
        sets.bar.color = color_palette_method,
        sets.bar.show.numbers = FALSE,
        sets.bar.x.label = "Set Size",
        intersection.matrix.color = "#21bda3",
        specific = TRUE
      )
    }

    venn_plots[[as.character(threshold)]] <- venn_plot
  }

  return(venn_plots)
}

create_weight_plots <- function(
    network_table,
    intersection_results,
    thresholds = 500,
    point_size = 1.5,
    color_palette = c("grey70", "#ac0000")) {
  weight_plots <- list()

  log_message("Creating weight plots for HARNexus...")
  network_table_infercsn <- network_table[["HARNexus"]]
  colnames(network_table_infercsn) <- c("TF", "gene", "weight")
  hic_intersect_genes <- intersection_results$intersect_genes_list$HARNexus[[as.character(0)]]

  network_table_infercsn$abs_weight <- abs(network_table_infercsn$weight)

  network_table_infercsn$is_hic_gene <- network_table_infercsn$gene %in% hic_intersect_genes

  network_table_infercsn <- network_table_infercsn[order(-network_table_infercsn$abs_weight), ]

  network_table_infercsn$rank <- 1:nrow(network_table_infercsn)

  weight_plots[["HARNexus"]] <- ggplot(
    network_table_infercsn,
    aes(x = rank, y = abs_weight, color = is_hic_gene)
  ) +
    geom_point(size = point_size) +
    scale_color_manual(
      values = color_palette,
      labels = c("Other edges", "Hi-C intersection edges")
    ) +
    labs(
      title = "HARNexus",
      x = "",
      y = "Absolute weight",
      color = "Edge type"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 30, hjust = 1)
    )

  network_table_3 <- list(
    "GENIE3" = network_table[["GENIE3"]],
    "LEAP" = network_table[["LEAP"]],
    "PPCOR" = network_table[["PPCOR"]]
  )

  for (method_name in names(network_table_3)) {
    for (threshold in thresholds) {
      log_message(
        "Creating weight plots for ",
        method_name, " of threshold ",
        threshold, "..."
      )
      weight_table <- network_table_3[[method_name]]
      colnames(weight_table) <- c("TF", "gene", "weight")
      weight_table_subset <- weight_table[1:threshold, ]
      hic_intersect_genes <- intersection_results$intersect_genes_list[[method_name]][[as.character(threshold)]]

      weight_table_subset$abs_weight <- abs(weight_table_subset$weight)

      weight_table_subset$is_hic_gene <- weight_table_subset$gene %in% hic_intersect_genes

      weight_table_subset <- weight_table_subset[order(-weight_table_subset$abs_weight), ]

      weight_table_subset$rank <- 1:nrow(weight_table_subset)

      p <- ggplot(
        weight_table_subset,
        aes(x = rank, y = abs_weight, color = is_hic_gene)
      ) +
        geom_point(size = point_size) +
        scale_color_manual(
          values = color_palette,
          labels = c("Other edges", "Hi-C intersection edges")
        ) +
        labs(
          title = method_name,
          x = "",
          y = "Absolute weight",
          color = "Edge type"
        ) +
        theme_bw() +
        theme(
          legend.position = "bottom",
          axis.text.x = element_text(angle = 30, hjust = 1)
        )

      weight_plots[[paste(method_name, threshold, sep = "_")]] <- p
    }
  }

  return(weight_plots)
}
