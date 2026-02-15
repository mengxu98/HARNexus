coverage_ratio <- function(
    network_table,
    method_name,
    step_size = 500,
    only_100 = TRUE) {
  all_target_genes <- unique(network_table[[2]])
  total_possible_genes <- length(all_target_genes)

  log_message(
    "{.pkg {method_name}} - Total target genes: {.val {total_possible_genes}}"
  )

  max_threshold <- nrow(network_table)
  thresholds <- seq(step_size, max_threshold, by = step_size)
  if (max(thresholds) < max_threshold) {
    thresholds <- c(thresholds, max_threshold)
  }

  coverage_data <- data.frame(
    threshold = numeric(),
    genes_found = numeric(),
    percentage = numeric()
  )

  found_all_genes <- FALSE
  min_threshold_all_genes <- NA

  for (threshold in thresholds) {
    current_subset <- network_table[1:threshold, ]
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
      log_message(
        "{.pkg {method_name}} - Found all genes at threshold: {.val {threshold}}"
      )

      if (only_100) {
        log_message(
          "{.pkg {method_name}} - Stopping at 100% coverage (only_100 = TRUE)"
        )
        break
      }
    }

    log_message(
      "{.pkg {method_name}} - Edges: {.val {threshold}}, Genes found: {.val {current_genes}}/{.val {total_possible_genes}} ({.val {percentage * 100}%})"
    )
  }

  return(
    list(
      method = method_name,
      coverage_data = coverage_data,
      total_possible = total_possible_genes,
      min_threshold_all_genes = min_threshold_all_genes,
      found_all_genes = found_all_genes
    )
  )
}

coverage_trend_plot <- function(
    network_list,
    color_methods,
    step_size = 500,
    only_100 = TRUE,
    text_size = 3) {
  genie3_coverage <- coverage_ratio(
    network_list[["GENIE3"]],
    "GENIE3",
    step_size = step_size,
    only_100 = only_100
  )
  harnexus_coverage <- coverage_ratio(
    network_list[["HARNexus"]],
    "HARNexus",
    step_size = step_size,
    only_100 = only_100
  )
  leap_coverage <- coverage_ratio(
    network_list[["LEAP"]],
    "LEAP",
    step_size = step_size,
    only_100 = only_100
  )
  ppcor_coverage <- coverage_ratio(
    network_list[["PPCOR"]],
    "PPCOR",
    step_size = step_size,
    only_100 = only_100
  )

  plot_data <- rbind(
    cbind(Method = "GENIE3", genie3_coverage$coverage_data),
    cbind(Method = "HARNexus", harnexus_coverage$coverage_data),
    cbind(Method = "LEAP", leap_coverage$coverage_data),
    cbind(Method = "PPCOR", ppcor_coverage$coverage_data)
  )

  final_summary <- data.frame(
    Method = c("GENIE3", "HARNexus", "LEAP", "PPCOR"),
    Total_Edges = c(
      nrow(network_list[["GENIE3"]]),
      nrow(network_list[["HARNexus"]]),
      nrow(network_list[["LEAP"]]),
      nrow(network_list[["PPCOR"]])
    ),
    Total_Genes = c(
      genie3_coverage$total_possible,
      harnexus_coverage$total_possible,
      leap_coverage$total_possible,
      ppcor_coverage$total_possible
    ),
    min_edges_all_genes = c(
      genie3_coverage$min_threshold_all_genes,
      harnexus_coverage$min_threshold_all_genes,
      leap_coverage$min_threshold_all_genes,
      ppcor_coverage$min_threshold_all_genes
    ),
    Max_coverage = c(
      max(genie3_coverage$coverage_data$percentage),
      max(harnexus_coverage$coverage_data$percentage),
      max(leap_coverage$coverage_data$percentage),
      max(ppcor_coverage$coverage_data$percentage)
    )
  )

  final_summary$max_coverage_percent <- sprintf(
    "%.2f%%", final_summary$Max_coverage * 100
  )
  final_summary$edges_ratio <- sprintf(
    "%.2f%%",
    final_summary$min_edges_all_genes / final_summary$Total_Edges * 100
  )

  ggplot(
    plot_data,
    aes(x = threshold, y = percentage, color = Method)
  ) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 0.9, alpha = 0.6) +
    geom_vline(
      data = final_summary,
      aes(xintercept = min_edges_all_genes, color = Method),
      linetype = "dashed",
      alpha = 0.7,
      show.legend = FALSE
    ) +
    geom_text(
      data = final_summary,
      aes(
        x = min_edges_all_genes,
        y = 0.8,
        label = sprintf("%d edges", min_edges_all_genes),
        color = Method
      ),
      angle = 90,
      vjust = 1.2,
      size = text_size,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = color_methods
    ) +
    scale_y_continuous(
      labels = scales::percent,
      limits = c(0, 1),
      breaks = seq(0, 1, 0.1)
    ) +
    scale_x_continuous(
      breaks = seq(0, max(plot_data$threshold), by = 2000)
    ) +
    labs(
      x = "Number of edges",
      y = "Genes coverage percentage",
      color = "Method"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

reorder_data <- function(values) {
  names(values) <- methods_order
  values[methods_order]
}

create_split_data <- function(methods_order, values, break_points) {
  df <- data.frame(
    Method = factor(names(values), levels = methods_order),
    Value = unname(values)
  )

  if (is.null(break_points)) {
    df$type <- "single"
    return(df)
  }

  df$type <- ifelse(df$Value > break_points[2], "big", "small")

  df_big <- filter(df, type == "big")
  if (nrow(df_big) > 0) {
    df_small_placeholder <- df_big %>%
      mutate(type = "small", Value = break_points[1], is_placeholder = TRUE)
    df$is_placeholder <- FALSE
    df_result <- bind_rows(df_small_placeholder, df)
  } else {
    df$is_placeholder <- FALSE
    df_result <- df
  }

  return(df_result)
}

create_split_plot <- function(
    data,
    title,
    y_lab,
    fill_color,
    break_points = NULL,
    legend_position = "right") {
  p <- ggplot(data, aes(x = Method, y = Value, fill = Method))
  if (is.null(break_points)) {
    p <- p +
      geom_rect(
        aes(
          xmin = as.numeric(Method) - 0.3,
          xmax = as.numeric(Method) + 0.3,
          ymin = 0,
          ymax = Value
        ),
        colour = "black",
        linewidth = 0.5
      ) +
      geom_text(
        aes(label = format(Value, big.mark = ",")),
        vjust = -0.5, size = 3
      )
  } else {
    if (!"is_placeholder" %in% names(data)) {
      data$is_placeholder <- FALSE
    }
    data$ymin <- ifelse(data$type == "small", 0, break_points[2])
    data$show_label <- !data$is_placeholder
    p <- p +
      geom_rect(
        aes(
          xmin = as.numeric(Method) - 0.3,
          xmax = as.numeric(Method) + 0.3,
          ymin = ymin,
          ymax = Value
        ),
        colour = "black",
        linewidth = 0.5
      ) +
      geom_text(
        aes(label = ifelse(show_label, format(Value, big.mark = ","), "")),
        vjust = -0.5, size = 3
      ) +
      facet_grid(type ~ ., scales = "free_y")
  }
  p <- p +
    scale_fill_manual(values = fill_color) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = legend_position,
      axis.text.x = element_text(angle = 30, hjust = 1)
    ) +
    labs(title = title, y = y_lab, x = "") +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.20)),
      labels = function(x) {
        format(x, scientific = FALSE, big.mark = ",", digits = 1)
      }
    )

  return(p)
}
