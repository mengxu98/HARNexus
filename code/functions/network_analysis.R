calculate_jaccard <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  return(intersection / union)
}

calculate_jaccard_intersect <- function(
    stage_pairs_list,
    stages,
    intersect_tfs,
    intersect_targets) {
  n_stages <- length(stages)
  jaccard_matrix <- matrix(0, nrow = n_stages, ncol = n_stages)
  rownames(jaccard_matrix) <- stages
  colnames(jaccard_matrix) <- stages

  processed_pairs <- lapply(stage_pairs_list, function(pairs) {
    split_pairs <- strsplit(pairs, "-")
    valid_pairs <- sapply(split_pairs, function(pair) {
      pair[1] %in% intersect_tfs && pair[2] %in% intersect_targets
    })
    pairs[valid_pairs]
  })

  for (i in 1:n_stages) {
    for (j in 1:n_stages) {
      jaccard_matrix[i, j] <- calculate_jaccard(
        processed_pairs[[stages[i]]],
        processed_pairs[[stages[j]]]
      )
    }
  }

  return(jaccard_matrix)
}


plot_heatmap <- function(
    jaccard_matrix,
    color_palette = c("white", "#961c45")) {
  col_fun <- circlize::colorRamp2(
    c(0, 1),
    color_palette
  )
  ht <- ComplexHeatmap::Heatmap(
    jaccard_matrix,
    name = "Jaccard\nsimilarity",
    col = col_fun,
    border = "gray10",
    rect_gp = grid::gpar(col = "gray50"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_names_rot = 0,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.2f", jaccard_matrix[i, j]), x, y)
    },
    column_title = "Developmental stages",
    row_title = "Developmental stages",
    heatmap_legend_param = list(
      title = "Jaccard\nsimilarity",
      at = c(0, 0.5, 1),
      labels = c("0", "0.5", "1")
    ),
    width = 4.5,
    height = 4.5
  )
  draw(ht)
}

upset_plot <- function(data) {
  data_upset <- as.data.frame(
    ComplexHeatmap::list_to_matrix(data)
  )
  data_upset$intersect <- rowSums(data_upset)
  highlight_color <- "#3333ff"
  queries_list <- list(
    list(
      query = elements,
      params = list("S3"),
      color = highlight_color,
      active = TRUE,
      query.name = "Diff"
    ),
    list(
      query = intersects,
      params = list("S4"),
      color = highlight_color,
      active = TRUE,
      query.name = "Diff"
    ),
    list(
      query = intersects,
      params = list("S5"),
      color = highlight_color,
      active = TRUE,
      query.name = "Diff"
    ),
    list(
      query = intersects,
      params = list("S6"),
      color = highlight_color,
      active = TRUE,
      query.name = "Diff"
    ),
    list(
      query = intersects,
      params = list("S7"),
      color = highlight_color,
      active = TRUE,
      query.name = "Diff"
    ),
    list(
      query = intersects,
      params = list("S8"),
      color = highlight_color,
      active = TRUE,
      query.name = "Diff"
    ),
    list(
      query = intersects,
      params = list("S9"),
      color = highlight_color,
      active = TRUE,
      query.name = "Diff"
    )
  )

  UpSetR::upset(
    data_upset,
    order.by = "freq",
    queries = queries_list,
    nsets = length(target_lists),
    nintersects = 30,
    sets.bar.color = "black",
    sets = NULL,
    matrix.color = "gray50",
    main.bar.color = "black",
    mainbar.y.label = "Intersection size",
    sets.x.label = "Set size",
    point.size = 3.5, line.size = 1,
    mb.ratio = c(0.6, 0.4),
    att.pos = NULL,
    show.numbers = "yes", number.angles = 0,
    group.by = "degree",
    cutoff = NULL,
    shade.color = "gray70",
    shade.alpha = 0.25,
    matrix.dot.alpha = 0.5,
    empty.intersections = NULL,
    color.pal = 2,
    attribute.plots = NULL, scale.intersections = "identity",
    scale.sets = "identity", text.scale = 1, set_size.angles = 0,
    set_size.show = FALSE, set_size.numbers_size = NULL,
    set_size.scale_max = NULL
  )
}

calculate_enrichment <- function(
    background_genes,
    target_genes,
    disease_genes) {
  in_disease <- sum(target_genes %in% disease_genes)
  not_in_disease <- sum(!(target_genes %in% disease_genes))
  in_background <- sum(background_genes %in% disease_genes)
  not_in_background <- sum(!(background_genes %in% disease_genes))

  contingency_table <- matrix(
    c(
      in_disease,
      not_in_disease,
      in_background - in_disease,
      not_in_background - not_in_disease
    ),
    nrow = 2,
    byrow = TRUE
  )

  fisher_test <- fisher.test(contingency_table)

  observed_ratio <- in_disease / length(target_genes)
  expected_ratio <- in_background / length(background_genes)
  fold_enrichment <- observed_ratio / expected_ratio

  return(
    list(
      p_value = fisher_test$p.value,
      fold_enrichment = fold_enrichment,
      observed = in_disease,
      expected = length(target_genes) * expected_ratio
    )
  )
}

create_layout_with_umap <- function(g, weights = NULL) {
  adj_matrix <- as.matrix(igraph::as_adjacency_matrix(g, attr = "Weight"))

  set.seed(42)
  umap_coords <- uwot::umap(
    adj_matrix,
    n_components = 2,
    n_neighbors = 15,
    min_dist = 0.1
  )

  layout <- data.frame(
    x = umap_coords[, 1],
    y = umap_coords[, 2]
  )

  layout$x <- scale(layout$x, center = TRUE, scale = TRUE)
  layout$y <- scale(layout$y, center = TRUE, scale = TRUE)

  return(as.matrix(layout))
}

plot_networks <- function(
    stage_data,
    stage_name,
    tfs,
    targets,
    max_edges_per_type = 500,
    split = FALSE,
    layout_type = "fr",
    label_size = 3.5,
    show_tf_labels = TRUE,
    show_gene_labels = TRUE,
    node_colors = c(TF = "#d41313", Gene = "#065e91"),
    node_size = c(TF = 3, Gene = 3),
    node_color_celltypes = FALSE) {
  cell_type_colors <- c(
    "Astro" = "#005ea3",
    "Endo" = "#24B700",
    "Micro" = "#00C1AB",
    "OPC" = "#00ACFC",
    "ExN" = "#f0749d",
    "InN" = "#c21b90",
    "NPC" = "#e29828",
    "Olig" = "#5865d3",
    "Perc" = "#c08f09"
  )
  edges <- stage_data[
    stage_data$TF %in% tfs &
      stage_data$Target %in% targets,
    c("TF", "Target", "Weight", "CellType")
  ]

  if (nrow(edges) == 0) {
    return(NULL)
  }

  edges <- do.call(rbind, lapply(unique(edges$CellType), function(ct) {
    ct_edges <- edges[edges$CellType == ct, ]
    ct_edges <- ct_edges[order(abs(ct_edges$Weight), decreasing = TRUE), ]
    head(ct_edges, max_edges_per_type)
  }))

  g_full <- graph_from_data_frame(edges, directed = TRUE)
  V(g_full)$type <- ifelse(V(g_full)$name %in% tfs, "TF", "Gene")
  V(g_full)$celltype <- edges$CellType[match(V(g_full)$name, edges$Target)]
  V(g_full)$celltype[V(g_full)$type == "TF"] <- "TF"

  V(g_full)$importance <- sapply(
    V(g_full)$name, function(node) {
      if (node %in% tfs) {
        out_edges <- incident(g_full, node, mode = "out")
        sum(abs(E(g_full)[out_edges]$Weight))
      } else {
        in_edges <- incident(g_full, node, mode = "in")
        sum(abs(E(g_full)[in_edges]$Weight))
      }
    }
  )

  V(g_full)$size <- ifelse(
    V(g_full)$type == "TF",
    node_size["TF"],
    node_size["Gene"]
  )

  cell_type_order <- unique(edges$CellType)

  missing_types <- setdiff(cell_type_order, names(cell_type_colors))
  if (length(missing_types) > 0) {
    additional_colors <- scales::hue_pal()(length(missing_types))
    names(additional_colors) <- missing_types
    cell_type_colors <- c(cell_type_colors, additional_colors)
  }

  cell_type_colors <- cell_type_colors[cell_type_order]

  if (split) {
    plots <- lapply(names(cell_type_colors), function(ct) {
      ct_edges <- edges[edges$CellType == ct, ]
      g_ct <- graph_from_data_frame(ct_edges, directed = TRUE)

      V(g_ct)$type <- V(g_full)$type[match(V(g_ct)$name, V(g_full)$name)]
      V(g_ct)$importance <- V(g_full)$importance[match(V(g_ct)$name, V(g_full)$name)]
      V(g_ct)$size <- ifelse(
        V(g_ct)$type == "TF",
        node_size["TF"],
        node_size["Gene"]
      )

      set.seed(42)
      if (layout_type == "umap") {
        layout <- create_layout_with_umap(g_ct)
        layout <- as.data.frame(layout)
        colnames(layout) <- c("x", "y")
      } else {
        layout <- create_layout(
          g_ct,
          layout = layout_type,
          weights = abs(E(g_ct)$Weight),
          maxiter = 1000
        )
      }

      layout$x <- layout$x * 1.5
      layout$y <- layout$y * 1.5

      p <- ggraph(g_ct, layout = layout) +
        geom_edge_link(
          aes(
            edge_width = abs(Weight) * 100,
            edge_colour = Weight
          ),
          end_cap = circle(1, "mm")
        ) +
        scale_edge_width_continuous(
          range = c(0.05, 0.3),
          name = "Weight"
        ) +
        geom_node_point(
          aes(
            size = size,
            shape = type,
            color = type
          )
        ) +
        scale_shape_manual(
          values = c("TF" = "diamond", "Gene" = "circle"),
          name = "Type"
        ) +
        scale_color_manual(
          values = node_colors,
          name = "Type"
        ) +
        scale_size_identity()

      if (show_tf_labels || show_gene_labels) {
        p <- p + geom_node_text(
          data = function(x) {
            show_nodes <- (show_tf_labels & x$type == "TF") |
              (show_gene_labels & x$type == "Gene" &
                x$importance > quantile(x$importance[x$type == "Gene"], 0.95))
            x[show_nodes, ]
          },
          aes(label = name),
          repel = TRUE,
          size = label_size,
          fontface = "bold",
          max.overlaps = 20,
          bg.colour = "white",
          bg.r = 0.15
        )
      }

      p <- p + theme_graph(
        base_family = "sans",
        background = "white"
      ) +
        theme(
          legend.position = "right",
          legend.box = "vertical"
        ) +
        labs(
          title = sprintf("%s - %s", stage_name, ct)
        )

      return(p)
    })

    combined_plot <- wrap_plots(
      plots,
      ncol = 3,
      guides = "collect"
    )

    return(combined_plot)
  } else {
    set.seed(42)
    if (layout_type == "umap") {
      layout <- create_layout_with_umap(g_full)
      layout <- as.data.frame(layout)
      colnames(layout) <- c("x", "y")
    } else {
      layout <- create_layout(
        g_full,
        layout = layout_type,
        weights = abs(E(g_full)$Weight),
        maxiter = 1000
      )
    }

    layout$x <- layout$x * 1.5
    layout$y <- layout$y * 1.5

    p <- ggraph(g_full, layout = layout) +
      geom_edge_link(
        aes(
          edge_width = abs(Weight),
          edge_colour = if (!node_color_celltypes) CellType else "gray50"
        ),
        end_cap = circle(1, "mm")
      ) +
      scale_edge_width_continuous(
        range = c(0.01, 1),
        name = "Weight"
      )

    if (!node_color_celltypes) {
      p <- p + scale_edge_colour_manual(
        values = cell_type_colors,
        name = "Cell Type"
      )
    } else {
      p <- p + scale_edge_colour_manual(
        values = "gray50",
        guide = "none"
      )
    }

    p <- p + geom_node_point(
      aes(
        size = size,
        shape = type,
        color = if (!node_color_celltypes) type else celltype
      )
    ) +
      scale_shape_manual(
        values = c("TF" = "diamond", "Gene" = "circle"),
        name = "Type"
      )

    if (!node_color_celltypes) {
      p <- p + scale_color_manual(
        values = node_colors,
        name = "Type"
      )
    } else {
      node_colors_with_celltypes <- c("TF" = "black", cell_type_colors)
      p <- p + scale_color_manual(
        values = node_colors_with_celltypes,
        name = "Type"
      )
    }

    p <- p + scale_size_identity()

    if (show_tf_labels || show_gene_labels) {
      p <- p + geom_node_text(
        data = function(x) {
          show_nodes <- (show_tf_labels & x$type == "TF") |
            (show_gene_labels & x$type == "Gene" &
              x$importance > quantile(x$importance[x$type == "Gene"], 0.95))
          x[show_nodes, ]
        },
        aes(label = name),
        repel = TRUE,
        size = label_size,
        fontface = "bold",
        max.overlaps = 20,
        bg.colour = "white",
        bg.r = 0.15
      )
    }

    p <- p + theme_graph(
      base_family = "sans",
      background = "white"
    ) +
      theme(
        legend.position = "right",
        legend.box = "vertical",
        plot.margin = margin(10, 10, 10, 10)
      ) +
      labs(
        title = paste("PFC", stage_name, "network")
      )

    return(p)
  }
}

plot_wordcloud <- function(
    network_data,
    type = c("TF", "Target"),
    genes = NULL,
    min_freq = 3,
    max_words = 100,
    random_order = FALSE,
    colors = c(
      "#1B9E77",
      "#105e1d",
      "#7570B3",
      "#E7298A",
      "#db6f31",
      "#305df1"
    ),
    shape = "square") {
  type <- match.arg(type)

  gene_stats <- network_data %>%
    filter(if (type == "TF") TF %in% genes else Target %in% genes) %>%
    group_by(if (type == "TF") TF else Target) %>%
    summarise(
      weight_sum = sum(abs(Weight)),
      frequency = n(),
      mean_weight = mean(abs(Weight)),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_weight))

  colnames(gene_stats)[1] <- "gene"

  word_freq <- data.frame(
    word = gene_stats$gene,
    freq = scale(gene_stats$mean_weight)[, 1] * 10 + 20
  )

  set.seed(42)
  p <- ggplot(
    word_freq,
    aes(
      label = word,
      size = freq,
      color = freq
    )
  ) +
    ggwordcloud::geom_text_wordcloud(
      shape = shape,
      rm_outside = TRUE,
      grid_size = 1,
      max_steps = 1000,
      eccentricity = 1
    ) +
    scale_size_area(max_size = 24) +
    scale_color_gradientn(colors = colors) +
    theme_minimal() +
    theme(
      plot.background = element_rect(
        fill = "white",
        color = "black",
        linewidth = 1
      ),
      plot.margin = margin(10, 10, 10, 10)
    )

  return(
    list(
      plot = p,
      stats = gene_stats
    )
  )
}

plot_importance <- function(
    network_data,
    type = c("TF", "Target"),
    genes = NULL,
    top_n = 20,
    color_gradient = c("#7b9ecc", "#1d559e")) {
  type <- match.arg(type)

  gene_stats <- network_data %>%
    filter(if (type == "TF") {
      .data$TF %in%
        genes
    } else {
      .data$Target %in% genes
    }) %>%
    group_by(gene = if (type == "TF") .data$TF else .data$Target) %>%
    summarise(
      weight_sum = sum(abs(.data$Weight)),
      frequency = n(),
      mean_weight = mean(abs(.data$Weight)),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_weight))

  top_genes <- head(gene_stats, top_n)

  p <- ggplot(
    top_genes,
    aes(
      x = reorder(gene, mean_weight),
      y = mean_weight
    )
  ) +
    geom_bar(
      stat = "identity",
      aes(fill = mean_weight),
      width = 0.7
    ) +
    scale_fill_gradient(
      low = color_gradient[1],
      high = color_gradient[2],
      name = "Mean absolute weight"
    ) +
    coord_flip() +
    theme_bw() +
    theme(
      axis.text.x = element_text(
        angle = 30,
        hjust = 0.5,
        vjust = 0.5
      )
    ) +
    labs(
      x = if (type == "TF") "Transcription factor" else "Target gene",
      y = "Mean absolute weight"
    )

  return(
    list(
      plot = p,
      stats = gene_stats
    )
  )
}
