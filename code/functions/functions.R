plot_network_distribution <- function(
    network_table,
    compare_random = TRUE) {
  g_orig <- igraph::graph_from_data_frame(
    network_table,
    directed = TRUE
  )
  degrees_orig <- igraph::degree(g_orig, mode = "total")

  degree_freq <- table(degrees_orig)
  df_orig <- data.frame(
    k = as.numeric(names(degree_freq)),
    P_k = as.numeric(degree_freq) / sum(degree_freq)
  )

  dist_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.major = element_line(color = "grey95", size = 0.2),
      panel.grid.minor = element_blank(),
      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
      panel.border = element_rect(color = "black", fill = NA, size = 0.3),
      aspect.ratio = 1
    )

  p1 <- ggplot2::ggplot(df_orig, aes(x = log(k), y = log(P_k))) +
    geom_point(size = 1) +
    geom_smooth(
      method = "lm",
      se = FALSE,
      color = "grey50",
      size = 0.5,
      formula = y ~ x
    ) +
    dist_theme +
    labs(
      x = "log k",
      y = "log P(k)",
      title = "Degree distribution",
      tag = "a"
    )

  if (nrow(df_orig) > 1) {
    model <- lm(log(P_k) ~ log(k), data = df_orig)
    r2 <- summary(model)$r.squared
    p1 <- p1 + annotate("text",
      x = max(log(df_orig$k)) - 0.1,
      y = max(log(df_orig$P_k)) - 0.1,
      label = sprintf("R² = %.3f", r2),
      size = 3,
      color = "steelblue",
      hjust = 1
    )
  }

  if (compare_random) {
    n_edges <- igraph::ecount(g_orig)
    all_nodes <- unique(c(network_table$regulator, network_table$target))

    random_edges <- data.frame(
      regulator = sample(all_nodes, n_edges, replace = TRUE),
      target = sample(all_nodes, n_edges, replace = TRUE),
      weight = network_table$weight
    )

    g_random <- igraph::graph_from_data_frame(
      random_edges,
      directed = TRUE
    )
    degrees_random <- igraph::degree(g_random, mode = "total")

    degree_freq_random <- table(degrees_random)
    df_random <- data.frame(
      k = as.numeric(names(degree_freq_random)),
      P_k = as.numeric(degree_freq_random) / sum(degree_freq_random)
    )

    p2 <- ggplot2::ggplot(df_random, aes(x = log(k), y = log(P_k))) +
      geom_point(size = 1) +
      geom_smooth(
        method = "lm",
        se = FALSE,
        color = "grey50",
        size = 0.3,
        formula = y ~ x
      ) +
      dist_theme +
      labs(
        x = "log k",
        y = "log P(k)",
        title = "Degree distribution\nof randomized network",
        tag = "b"
      )

    if (nrow(df_random) > 1) {
      model_random <- lm(log(P_k) ~ log(k), data = df_random)
      r2_random <- summary(model_random)$r.squared
      p2 <- p2 + annotate("text",
        x = max(log(df_random$k)) - 0.1,
        y = max(log(df_random$P_k)) - 0.1,
        label = sprintf("R² = %.3f", r2_random),
        size = 3,
        color = "steelblue",
        hjust = 1
      )
    }

    return(list(original = p1, random = p2))
  }

  return(p1)
}
