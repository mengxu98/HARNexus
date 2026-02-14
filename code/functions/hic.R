calculate_intersection_trend <- function(
    network_list,
    tfdb_data,
    hic_data,
    exclude_tfs_in_genes = TRUE,
    max_edges = 9000,
    step_size = 500) {
  results_df <- data.frame(
    threshold = numeric(),
    method = character(),
    total_genes_num = numeric(),
    intersect_genes_num = numeric(),
    ratio = numeric()
  )

  thresholds <- seq(step_size, max_edges, by = step_size)

  intersect_genes_list <- list()
  tfs_list <- list()
  target_genes_list <- list()

  for (threshold in thresholds) {
    for (method_name in names(network_list)) {
      network_table <- network_list[[method_name]]
      colnames(network_table) <- c("TF", "gene", "weight")
      if (exclude_tfs_in_genes) {
        tfs <- unique(network_table$TF)
        network_table <- network_table[!network_table$gene %in% tfs, ]
      }
      current_threshold <- min(threshold, nrow(network_table))
      network_table_sub <- network_table[1:current_threshold, ]


      network <- merge(network_table_sub, tfdb_data, by = "TF")
      network <- network[, c(4, 5, 1, 2, 3)]

      compar_hic <- merge(network, hic_data, by = "HAR")

      filtered_df <- compar_hic[compar_hic[["gene"]] == compar_hic$hic_gene, ]

      total_genes_num <- length(unique(compar_hic[["gene"]]))
      intersect_genes <- unique(filtered_df[["gene"]])
      intersect_genes_list[[method_name]][[as.character(threshold)]] <- as.character(
        intersect_genes
      )
      tfs_list[[method_name]][[as.character(threshold)]] <- as.character(
        unique(compar_hic[["TF"]])
      )
      target_genes_list[[method_name]][[as.character(threshold)]] <- as.character(
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
      ratio <- sprintf(
        "%.2f%%",
        ratio * 100
      )
      log_message(
        "{.val {method_name}}: edges: {.val {threshold}}, intersect: {.val {intersect_genes_num}/{total_genes_num}} ({.val {ratio}})"
      )
    }
  }

  return(
    list(
      results = results_df,
      tfs_list = tfs_list,
      target_genes_list = target_genes_list,
      intersect_genes_list = intersect_genes_list
    )
  )
}

intersection_ratio_trend_plot <- function(
    intersection_results,
    color_methods,
    reference_method = "HARNexus") {
  results_df <- intersection_results$results
  methods <- unique(results_df$method)
  max_threshold <- max(results_df$threshold)

  final_points <- results_df[results_df$threshold == max_threshold, ]

  final_points$label <- sprintf(
    "%s: %.1f%%",
    final_points$method,
    final_points$ratio * 100
  )

  ref_ratio <- final_points$ratio[final_points$method == reference_method]
  other_methods <- setdiff(methods, reference_method)

  bracket_df <- data.frame(
    method = character(),
    fold = numeric(),
    y_ref = numeric(),
    y_other = numeric(),
    fold_label = character(),
    label_y = numeric(),
    stringsAsFactors = FALSE
  )

  for (m in other_methods) {
    other_ratio <- final_points$ratio[final_points$method == m]
    if (other_ratio > 0) {
      fold <- ref_ratio / other_ratio
      bracket_df <- rbind(
        bracket_df,
        data.frame(
          method = m,
          fold = fold,
          y_ref = ref_ratio,
          y_other = other_ratio,
          fold_label = sprintf("%.1fx", fold),
          label_y = (ref_ratio + other_ratio) / 2,
          stringsAsFactors = FALSE
        )
      )
    }
  }

  bracket_df <- bracket_df[order(bracket_df$fold), ]

  x_gap <- max_threshold * 0.07
  bracket_df$x_bracket <- max_threshold + x_gap * seq_len(nrow(bracket_df))

  x_max_bracket <- max(bracket_df$x_bracket)

  p <- ggplot(
    results_df,
    aes(x = threshold, y = ratio, color = method)
  ) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 0.8) +
    geom_point(
      data = final_points,
      aes(x = threshold, y = ratio, color = method),
      size = 2.5,
      shape = 21,
      fill = "white",
      stroke = 1.5
    ) +
    ggrepel::geom_text_repel(
      data = final_points,
      aes(x = threshold, y = ratio, label = label, color = method),
      size = 3,
      bg.color = "white",
      bg.r = 0.2,
      min.segment.length = 0,
      segment.size = 0.5,
      segment.linetype = "dashed",
      box.padding = 0.5,
      point.padding = 0.3,
      nudge_y = 0.015,
      direction = "y",
      force = 5,
      force_pull = 0.5,
      max.overlaps = Inf,
      xlim = c(NA, max_threshold * 0.85),
      show.legend = FALSE
    )

  p <- p + annotate(
    "segment",
    x = max_threshold, xend = x_max_bracket,
    y = ref_ratio, yend = ref_ratio,
    color = color_methods[reference_method],
    linetype = "dashed", linewidth = 0.4
  )

  tick_w <- max_threshold * 0.015

  for (i in seq_len(nrow(bracket_df))) {
    x_b <- bracket_df$x_bracket[i]

    p <- p + annotate(
      "segment",
      x = max_threshold, xend = x_b,
      y = bracket_df$y_other[i], yend = bracket_df$y_other[i],
      color = color_methods[bracket_df$method[i]],
      linetype = "dashed", linewidth = 0.4
    )

    p <- p + annotate(
      "segment",
      x = x_b, xend = x_b,
      y = bracket_df$y_ref[i], yend = bracket_df$y_other[i],
      color = color_methods[bracket_df$method[i]],
      linewidth = 0.5
    )

    p <- p + annotate(
      "segment",
      x = x_b - tick_w, xend = x_b + tick_w,
      y = bracket_df$y_ref[i], yend = bracket_df$y_ref[i],
      color = color_methods[bracket_df$method[i]],
      linewidth = 0.5
    )
    p <- p + annotate(
      "segment",
      x = x_b - tick_w, xend = x_b + tick_w,
      y = bracket_df$y_other[i], yend = bracket_df$y_other[i],
      color = color_methods[bracket_df$method[i]],
      linewidth = 0.5
    )

    p <- p + annotate(
      "text",
      x = x_b + max_threshold * 0.025,
      y = bracket_df$label_y[i],
      label = bracket_df$fold_label[i],
      color = color_methods[bracket_df$method[i]],
      size = 2.8, fontface = "bold", hjust = 0.5, angle = 90
    )
  }

  p <- p +
    labs(
      x = "Number of edges",
      y = "Intersection ratio",
      color = "Method"
    ) +
    theme_bw() +
    scale_color_manual(
      values = color_methods[names(color_methods) %in% results_df$method]
    ) +
    scale_y_continuous(
      labels = scales::percent,
      limits = c(0, max(results_df$ratio) * 1.05),
      breaks = seq(0, max(results_df$ratio), 0.1)
    ) +
    scale_x_continuous(
      breaks = seq(0, max_threshold, by = 1000),
      limits = c(0, x_max_bracket + max_threshold * 0.03)
    ) +
    coord_cartesian(clip = "off") +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 15, 5, 5, "pt")
    )

  p
}

comparison_barplot <- function(
    intersection_results,
    thresholds = c(500, 1000, 2000),
    color_methods,
    text_size = 2.5) {
  thresholds <- c(0, thresholds)
  filtered_data <- intersection_results$results[
    intersection_results$results$threshold %in% thresholds,
  ]

  filtered_data$method <- factor(
    filtered_data$method,
    levels = c("GENIE3", "HARNexus", "LEAP", "PPCOR")
  )

  ggplot() +
    geom_bar(
      data = filtered_data,
      aes(x = method, y = ratio, color = method),
      fill = NA,
      stat = "identity",
      width = 0.8
    ) +
    scale_color_manual(
      values = color_methods,
      name = "Method"
    ) +
    labs(
      title = "",
      x = "Method",
      y = "Intersection ratio",
      color = "Method"
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
      size = text_size
    ) +
    coord_flip() +
    theme_bw() +
    theme(
      legend.position = "right",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_y_continuous(
      limits = c(0, max(filtered_data$ratio) * 1.05),
      labels = scales::percent
    )
}

extract_genes_list <- function(
    results,
    type = c("intersect", "all", "tfs", "setdiff"),
    threshold = 500) {
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
      unique(results$intersect_genes_list$HARNexus[[as.character(threshold)]])
    )
  } else if (type == "all") {
    genes_list[["GENIE3"]] <- as.character(
      unique(results$target_genes_list$GENIE3[[as.character(threshold)]])
    )
    genes_list[["LEAP"]] <- as.character(
      unique(results$target_genes_list$LEAP[[as.character(threshold)]])
    )
    genes_list[["PPCOR"]] <- as.character(
      unique(results$target_genes_list$PPCOR[[as.character(threshold)]])
    )
    genes_list[["HARNexus"]] <- as.character(
      unique(results$target_genes_list$HARNexus[[as.character(threshold)]])
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
      unique(results$tfs_list$HARNexus[[as.character(threshold)]])
    )
  } else if (type == "setdiff") {
    genes_genie3_all <- as.character(
      unique(results$target_genes_list$GENIE3[[as.character(threshold)]])
    )
    genes_genie3_intersect <- as.character(
      unique(results$intersect_genes_list$GENIE3[[as.character(threshold)]])
    )
    genes_list[["GENIE3"]] <- as.character(
      setdiff(genes_genie3_all, genes_genie3_intersect)
    )
    genes_leap_all <- as.character(
      unique(results$target_genes_list$LEAP[[as.character(threshold)]])
    )
    genes_leap_intersect <- as.character(
      unique(results$intersect_genes_list$LEAP[[as.character(threshold)]])
    )
    genes_list[["LEAP"]] <- as.character(
      setdiff(genes_leap_all, genes_leap_intersect)
    )
    genes_ppcor_all <- as.character(
      unique(results$target_genes_list$PPCOR[[as.character(threshold)]])
    )
    genes_ppcor_intersect <- as.character(
      unique(results$intersect_genes_list$PPCOR[[as.character(threshold)]])
    )
    genes_list[["PPCOR"]] <- as.character(
      setdiff(genes_ppcor_all, genes_ppcor_intersect)
    )
    genes_infercsn_all <- as.character(
      unique(results$target_genes_list$HARNexus[[as.character(threshold)]])
    )
    genes_infercsn_intersect <- as.character(
      unique(results$intersect_genes_list$HARNexus[[as.character(threshold)]])
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
    color_methods,
    type = c("intersect", "all", "tfs", "setdiff"),
    return_upset = FALSE) {
  venn_plots <- list()

  for (threshold in thresholds) {
    log_message(paste("Creating Venn diagram for threshold", threshold, "..."))

    method_genes <- extract_genes_list(
      results,
      type = type,
      threshold = threshold
    )
    method_genes <- list(
      HARNexus = method_genes$HARNexus,
      GENIE3 = method_genes$GENIE3,
      PPCOR = method_genes$PPCOR,
      LEAP = method_genes$LEAP
    )
    color_methods <- c(
      "HARNexus" = color_methods["HARNexus"],
      "GENIE3" = color_methods["GENIE3"],
      "PPCOR" = color_methods["PPCOR"],
      "LEAP" = color_methods["LEAP"]
    )

    venn_plot <- ggVennDiagram::ggVennDiagram(
      method_genes,
      category.names = names(method_genes),
      label_alpha = 0,
      label_color = "white",
      label_size = 5,
      set_color = color_methods,
      edge_lty = "solid",
      edge_size = 1.5
    ) +
      theme(
        legend.position = "none"
      )

    venn_plot_build <- ggplot2::ggplot_build(venn_plot)
    for (i in seq_along(venn_plot_build$data)) {
      layer_data <- venn_plot_build$data[[i]]
      if ("label" %in% names(layer_data) && "x" %in% names(layer_data)) {
        category_labels <- layer_data$label %in% names(method_genes)
        if (any(category_labels, na.rm = TRUE)) {
          x_center <- mean(range(layer_data$x, na.rm = TRUE))
          y_center <- mean(range(layer_data$y, na.rm = TRUE))
          for (j in which(category_labels)) {
            venn_plot_build$data[[i]]$x[j] <- layer_data$x[j] +
              (x_center - layer_data$x[j]) * 0.2
            if ("y" %in% names(layer_data)) {
              venn_plot_build$data[[i]]$y[j] <- layer_data$y[j] +
                (y_center - layer_data$y[j]) * 0.1
            }
          }
        }
      }
    }
    venn_plot <- ggplot2::ggplot_gtable(venn_plot_build)
    venn_plot <- patchwork::wrap_elements(venn_plot)

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
        sets.bar.color = color_methods,
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
    network_list,
    intersection_results,
    thresholds = 500,
    point_size = 1,
    color_palette = c("grey70", "#15559A"),
    color_methods,
    tfdb_data = NULL,
    hic_data = NULL,
    top_n_validated = 5) {
  weight_plots <- list()
  top_n_validated_sets <- list()
  validated_pairs_sets <- list()

  if (!is.null(tfdb_data) && !is.null(hic_data)) {
    log_message("Preparing data...")

    tf_har_gene_map <- merge(
      tfdb_data[, c("TF", "HAR")],
      hic_data[, c("HAR", "hic_gene")],
      by = "HAR",
      all = FALSE
    )
    valid_tf_gene_pairs <- unique(
      paste(tf_har_gene_map$TF, tf_har_gene_map$hic_gene, sep = "_")
    )

    for (method_name in names(network_list)) {
      log_message("Processing {.pkg {method_name}}...")
      network_table <- network_list[[method_name]]
      colnames(network_table) <- c("TF", "gene", "weight")
      network_table$abs_weight <- abs(network_table$weight)

      network_pairs <- paste(network_table$TF, network_table$gene, sep = "_")
      validated_mask <- network_pairs %in% valid_tf_gene_pairs
      if (any(validated_mask)) {
        validated_df <- network_table[validated_mask, c("TF", "gene", "abs_weight")]
        validated_df <- validated_df[order(-validated_df$abs_weight), ]
        pair_ids <- paste(validated_df$TF, validated_df$gene, sep = "_")
        validated_df <- validated_df[!duplicated(pair_ids), ]
        validated_df <- validated_df[order(-validated_df$abs_weight), ]
        validated_pairs_sets[[method_name]] <- paste(
          validated_df$TF, validated_df$gene,
          sep = "_"
        )
        top_n_validated_sets[[method_name]] <- head(
          validated_df[, c("TF", "gene", "abs_weight")], top_n_validated
        )
      } else {
        validated_pairs_sets[[method_name]] <- character(0)
        top_n_validated_sets[[method_name]] <- NULL
      }
    }
  }

  for (threshold in thresholds) {
    for (method_name in names(network_list)) {
      network_table <- network_list[[method_name]]
      colnames(network_table) <- c("TF", "gene", "weight")

      network_table_sub <- network_table[1:threshold, ]
      network_table_sub$abs_weight <- abs(network_table_sub$weight)

      if (method_name %in% names(validated_pairs_sets)) {
        network_table_sub$pair_id <- paste(
          network_table_sub$TF, network_table_sub$gene,
          sep = "_"
        )
        network_table_sub$is_validated <- network_table_sub$pair_id %in% validated_pairs_sets[[method_name]]
        network_table_sub$pair_id <- NULL
      } else {
        network_table_sub$is_validated <- FALSE
      }

      hic_intersect_genes <- intersection_results$intersect_genes_list[[method_name]][[as.character(threshold)]]
      network_table_sub$is_hic_gene <- network_table_sub$gene %in% hic_intersect_genes

      network_table_sub <- network_table_sub[order(-network_table_sub$abs_weight), ]
      network_table_sub$rank <- seq_len(nrow(network_table_sub))

      network_table_sub$point_size <- ifelse(
        network_table_sub$is_validated,
        point_size,
        point_size * 0.2
      )

      network_table_sub$is_top_validated <- FALSE
      validated_edges_in_plot <- network_table_sub[network_table_sub$is_validated, ]

      if (nrow(validated_edges_in_plot) > 0) {
        validated_edges_in_plot <- validated_edges_in_plot[order(-validated_edges_in_plot$abs_weight), ]
        pair_ids <- paste(validated_edges_in_plot$TF, validated_edges_in_plot$gene, sep = "_")
        validated_edges_in_plot <- validated_edges_in_plot[!duplicated(pair_ids), ]
        top_n_in_plot <- head(validated_edges_in_plot, top_n_validated)

        if (nrow(top_n_in_plot) > 0) {
          top_n_pair_ids <- paste(top_n_in_plot$TF, top_n_in_plot$gene, sep = "_")
          network_table_sub$pair_id <- paste(network_table_sub$TF, network_table_sub$gene, sep = "_")
          network_table_sub$is_top_validated <- network_table_sub$pair_id %in% top_n_pair_ids
          network_table_sub$pair_id <- NULL
        }
      }

      p <- ggplot(
        network_table_sub,
        aes(x = rank, y = abs_weight, color = is_validated)
      ) +
        geom_point(aes(size = point_size)) +
        scale_size_identity(guide = "none") +
        scale_color_manual(
          values = color_palette,
          labels = c("Non-validated edges", "Validated edges")
        ) +
        labs(
          title = method_name,
          x = "Edges",
          y = "",
          color = "Edge type"
        ) +
        theme_bw() +
        theme(
          plot.title = element_text(color = color_methods[[method_name]]),
          axis.text.x = element_text(angle = 30, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )

      if (any(network_table_sub$is_top_validated)) {
        top_validated_data <- network_table_sub[network_table_sub$is_top_validated, ]
        p <- p +
          geom_point(
            data = top_validated_data,
            aes(x = rank, y = abs_weight),
            color = "#15559A",
            size = 1,
            shape = 21,
            fill = "white",
            stroke = 1
          ) +
          ggrepel::geom_text_repel(
            data = top_validated_data,
            aes(x = rank, y = abs_weight, label = paste0(TF, "-", gene)),
            color = "#15559A",
            size = 2.5,
            min.segment.length = 0,
            segment.size = 0.3,
            segment.color = "#15559A",
            segment.linetype = "dashed",
            box.padding = 0.3,
            point.padding = 0.2,
            max.overlaps = 10,
            show.legend = FALSE
          )
      }

      weight_plots[[paste(method_name, threshold, sep = "_")]] <- p
    }
  }

  return(weight_plots)
}

create_ratio_miniplot <- function(
    intersection_results,
    threshold,
    method_name,
    color_methods,
    plot_type = c("bar", "pie"),
    text_size = 3,
    color_palette = c("grey50", "#15559A")) {
  names(color_palette) <- c("Non-validated", "Validated")
  plot_type <- match.arg(plot_type)

  filtered_data <- intersection_results$results[
    intersection_results$results$threshold == threshold &
      intersection_results$results$method == method_name,
  ]

  if (nrow(filtered_data) == 0) {
    return(NULL)
  }

  ratio <- filtered_data$ratio[1]
  validated_ratio <- ratio
  non_validated_ratio <- 1 - ratio

  if (plot_type == "bar") {
    plot_data <- data.frame(
      category = c("Validated", "Non-validated"),
      value = c(validated_ratio, non_validated_ratio)
    )

    p <- ggplot(
      plot_data, aes(x = 1, y = value, fill = category)
    ) +
      geom_bar(stat = "identity", width = 0.5) +
      scale_fill_manual(
        values = color_palette,
        guide = "none"
      ) +
      labs(x = "", y = "") +
      geom_text(
        aes(
          x = 1,
          y = validated_ratio / 2,
          label = sprintf("%.1f%%", validated_ratio * 100)
        ),
        color = "black",
        size = text_size
      ) +
      theme_bw() +
      theme(
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(2, 2, 2, 2, "pt"),
        panel.border = element_rect(color = "white", linewidth = 0.5)
      ) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_x_continuous(limits = c(0.5, 1.5))
  } else {
    plot_data <- data.frame(
      category = c("Validated", "Non-validated"),
      value = c(validated_ratio, non_validated_ratio)
    )

    plot_data$ymax <- cumsum(plot_data$value)
    plot_data$ymin <- c(0, head(plot_data$ymax, n = -1))
    plot_data$labelPosition <- (plot_data$ymax + plot_data$ymin) / 2

    p <- ggplot(
      plot_data,
      aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = category)
    ) +
      geom_rect() +
      coord_polar(theta = "y") +
      xlim(c(3, 4)) +
      scale_fill_manual(
        values = color_palette,
        guide = "none"
      ) +
      geom_text(
        aes(
          x = 3.5, y = labelPosition,
          label = sprintf("%.1f%%", value * 100)
        ),
        color = "white",
        size = text_size
      ) +
      theme_void() +
      theme(
        plot.margin = margin(2, 2, 2, 2, "pt"),
        plot.background = element_rect(fill = "white", color = NA)
      )
  }

  return(p)
}

go_enrichment <- function(
    gene_list,
    validated_genes = NULL,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    top_n = 10,
    ont = c("BP", "MF", "CC", "ALL")) {
  ont <- match.arg(ont)
  if (length(gene_list) == 0) {
    return(NULL)
  }

  tryCatch(
    {
      gene_ids <- bitr(
        gene_list,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db
      )

      if (nrow(gene_ids) == 0) {
        return(NULL)
      }

      go_result <- clusterProfiler::enrichGO(
        gene = gene_ids$ENTREZID,
        OrgDb = org.Hs.eg.db,
        ont = ont,
        pAdjustMethod = "BH",
        pvalueCutoff = pvalueCutoff,
        qvalueCutoff = qvalueCutoff,
        readable = TRUE
      )

      if (is.null(go_result) || nrow(go_result) == 0) {
        return(NULL)
      }

      go_df <- as.data.frame(go_result)
      go_df <- go_df[order(go_df$p.adjust), ]
      go_df <- head(go_df, top_n)

      return(go_df)
    },
    error = function(e) {
      log_message("GO enrichment failed: {.val {e$message}}")
      return(NULL)
    }
  )
}

enrichment_heatmap <- function(
    go_result,
    method_name,
    color_methods,
    validated_genes = NULL,
    other_term_color = "grey70",
    top_n = 10) {
  if (is.null(go_result) || nrow(go_result) == 0) {
    p <- ggplot() +
      annotate(
        "text",
        x = 0.5, y = 0.5, label = "No significant GO terms",
        size = 3,
        color = other_term_color
      ) +
      theme_void()
    return(p)
  }

  go_df <- go_result[order(go_result$p.adjust), ]
  go_df <- head(go_df, top_n)

  go_df$Description <- sapply(go_df$Description, thisutils::capitalize)

  all_genes <- unique(unlist(strsplit(go_df$geneID, "/")))
  if (!is.null(validated_genes) && length(validated_genes) > 0) {
    validated_genes_set <- validated_genes
  } else {
    validated_genes_set <- character(0)
  }

  heatmap_data <- data.frame()
  for (i in seq_len(nrow(go_df))) {
    pathway_genes <- strsplit(go_df$geneID[i], "/")[[1]]
    for (gene in pathway_genes) {
      if (gene %in% all_genes) {
        is_validated <- gene %in% validated_genes_set
        heatmap_data <- rbind(
          heatmap_data,
          data.frame(
            Pathway = go_df$Description[i],
            Gene = gene,
            Value = 1,
            is_validated = is_validated,
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }

  if (nrow(heatmap_data) == 0) {
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No genes in pathways") +
      theme_void()
    return(p)
  }

  pathway_order <- go_df$Description[order(go_df$p.adjust)]
  gene_counts <- table(heatmap_data$Gene)
  gene_order <- names(sort(gene_counts, decreasing = TRUE))

  heatmap_data$Pathway <- factor(heatmap_data$Pathway, levels = pathway_order)
  heatmap_data$Gene <- factor(heatmap_data$Gene, levels = rev(gene_order))

  color_val <- color_methods[[method_name]]
  p <- ggplot(
    heatmap_data,
    aes(x = Pathway, y = Gene)
  ) +
    geom_tile(
      aes(fill = ifelse(is_validated, "validated", "other")),
      color = "white",
      linewidth = 0.5
    ) +
    scale_fill_manual(
      values = c("validated" = color_val, "other" = other_term_color),
      guide = "none"
    ) +
    labs(
      x = "",
      y = ""
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(
        angle = 30,
        hjust = 1,
        size = 9,
        color = "black"
      ),
      # axis.text.y = element_text(
      #   size = 5
      # ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 5, 5, 5, "pt")
    )

  if (length(validated_genes_set) > 0) {
    gene_levels <- levels(heatmap_data$Gene)
    gene_colors <- ifelse(gene_levels %in% validated_genes_set, color_val, "black")

    p <- p + theme(
      axis.text.y = element_text(
        size = 6,
        color = gene_colors
      )
    )
  }

  return(p)
}

weight_ratio_plots <- function(
    network_list,
    intersection_results,
    thresholds = 500,
    point_size = 1,
    color_palette = c("grey70", "black"),
    color_methods,
    tfdb_data = NULL,
    hic_data = NULL,
    ratio_plot_type = c("pie", "bar"),
    top_n_validated = 5,
    inset_position = list(
      left = 0.5, bottom = 0.5, right = 0.98, top = 0.98
    ),
    pie_text_size = 4,
    ncol = 4,
    perform_go_enrichment = TRUE,
    pvalue_cutoff = 0.05,
    qvalue_cutoff = 0.2,
    go_top_n = 10,
    heatmap_position = c("bottom", "right", "top", "left"),
    brain_pattern = NULL,
    go_ontology = c("BP", "MF", "CC", "ALL")) {
  ratio_plot_type <- match.arg(ratio_plot_type)
  heatmap_position <- match.arg(heatmap_position)
  go_ontology <- match.arg(go_ontology)

  weight_plots_list <- create_weight_plots(
    network_list = network_list,
    intersection_results = intersection_results,
    thresholds = thresholds,
    point_size = point_size,
    color_palette = color_palette,
    color_methods = color_methods,
    tfdb_data = tfdb_data,
    hic_data = hic_data,
    top_n_validated = top_n_validated
  )

  method_names <- names(network_list)
  ratio_plots_list <- list()
  for (method_name in method_names) {
    ratio_plots_list[[method_name]] <- create_ratio_miniplot(
      intersection_results = intersection_results,
      threshold = thresholds,
      method_name = method_name,
      color_methods = color_methods,
      plot_type = ratio_plot_type,
      text_size = pie_text_size,
      color_palette = color_palette
    )
  }

  weight_plots_with_ratio <- list()
  for (method_name in method_names) {
    plot_key <- paste0(method_name, "_", thresholds)
    if (plot_key %in% names(weight_plots_list)) {
      main_plot <- weight_plots_list[[plot_key]]
      if (method_name %in% names(ratio_plots_list)) {
        ratio_plot <- ratio_plots_list[[method_name]]
        weight_plots_with_ratio[[method_name]] <- main_plot +
          patchwork::inset_element(
            ratio_plot,
            left = inset_position$left,
            bottom = inset_position$bottom,
            right = inset_position$right,
            top = inset_position$top,
            align_to = "plot"
          )
      } else {
        weight_plots_with_ratio[[method_name]] <- main_plot
      }
    }
  }

  weight_plots_with_ratio <- weight_plots_with_ratio[method_names]

  if (length(weight_plots_with_ratio) > 0) {
    first_method <- method_names[1]
    weight_plots_with_ratio[[first_method]] <- weight_plots_with_ratio[[first_method]] +
      labs(y = "Absolute weight")
  }

  if (perform_go_enrichment) {
    log_message(paste0("Performing GO ", go_ontology, " enrichment analysis..."))

    if (!is.null(tfdb_data) && !is.null(hic_data)) {
      tf_har_gene_map <- merge(
        tfdb_data[, c("TF", "HAR")],
        hic_data[, c("HAR", "hic_gene")],
        by = "HAR",
        all = FALSE
      )
      valid_tf_gene_pairs <- unique(
        paste(tf_har_gene_map$TF, tf_har_gene_map$hic_gene, sep = "_")
      )
    } else {
      valid_tf_gene_pairs <- character(0)
    }

    go_plots_list <- list()
    for (method_name in method_names) {
      log_message("Processing GO enrichment for {.pkg {method_name}}...")
      genes_list <- extract_genes_list(
        intersection_results,
        type = "all",
        threshold = thresholds
      )
      method_genes <- genes_list[[method_name]]

      validated_genes <- NULL
      if (length(valid_tf_gene_pairs) > 0) {
        network_table <- network_list[[method_name]]
        colnames(network_table) <- c("TF", "gene", "weight")

        current_threshold <- min(thresholds, nrow(network_table))
        network_table_sub <- network_table[1:current_threshold, ]

        network_pairs <- paste(
          network_table_sub$TF, network_table_sub$gene,
          sep = "_"
        )
        validated_mask <- network_pairs %in% valid_tf_gene_pairs
        if (any(validated_mask)) {
          validated_genes <- unique(
            network_table_sub$gene[validated_mask]
          )
        }
      }

      go_result <- go_enrichment(
        gene_list = method_genes,
        validated_genes = validated_genes,
        pvalueCutoff = pvalue_cutoff,
        qvalueCutoff = qvalue_cutoff,
        top_n = go_top_n,
        ont = go_ontology
      )

      if (!is.null(brain_pattern) && !is.null(go_result) && nrow(go_result) > 0) {
        go_result <- go_result[
          grepl(brain_pattern, go_result$Description, ignore.case = TRUE), ,
          drop = FALSE
        ]
        if (nrow(go_result) == 0) {
          go_result <- NULL
        }
      }

      go_plot <- enrichment_heatmap(
        go_result = go_result,
        method_name = method_name,
        color_methods = color_methods,
        validated_genes = validated_genes,
        top_n = go_top_n
      )
      go_plots_list[[method_name]] <- go_plot
    }

    combined_plots_list <- list()
    for (method_name in method_names) {
      weight_plot <- weight_plots_with_ratio[[method_name]]
      go_plot <- go_plots_list[[method_name]]

      if (heatmap_position == "right") {
        combined_plots_list[[method_name]] <- weight_plot | go_plot +
          plot_layout(widths = c(0.6, 0.4))
      } else if (heatmap_position == "bottom") {
        combined_plots_list[[method_name]] <- weight_plot / go_plot +
          plot_layout(heights = c(0.2, 0.8))
      } else if (heatmap_position == "top") {
        combined_plots_list[[method_name]] <- go_plot / weight_plot +
          plot_layout(heights = c(0.8, 0.2))
      } else if (heatmap_position == "left") {
        combined_plots_list[[method_name]] <- go_plot | weight_plot +
          plot_layout(widths = c(0.4, 0.6))
      }
    }

    combined_plot <- patchwork::wrap_plots(combined_plots_list) +
      patchwork::plot_layout(
        ncol = ncol,
        guides = "collect"
      ) &
      theme(
        legend.position = "bottom"
      )
  } else {
    combined_plot <- patchwork::wrap_plots(weight_plots_with_ratio) +
      patchwork::plot_layout(
        ncol = ncol,
        guides = "collect"
      ) &
      theme(
        legend.position = "bottom"
      )
  }

  return(combined_plot)
}
