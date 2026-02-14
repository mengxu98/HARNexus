source("code/functions/prepare_env.R")

sample_pairs <- list(
  list(human = "h4", chimp = "c4"),
  list(human = "h3", chimp = "c2"),
  list(human = "h1", chimp = "c1")
)

human_color <- "#3271AE"
chimp_color <- "#D11A2D"
both_color <- "black"

test_show <- "Regulatory rewiring"
test_show <- NULL
top_n <- 0

trajectory_colors <- c(
  "Regulatory rewiring" = "#006D87",
  "Regulatory innovation" = "#F9BD10",
  "Cis-regulatory activation" = "#0AA344",
  "Conserved" = "#5E7987",
  "TF" = "#000000"
)

selected_celltypes <- c(
  "Astrocytes", "Excitatory neurons", "Inhibitory neurons",
  "Oligodendrocyte progenitor cells", "Microglia"
)

for (pair in sample_pairs) {
  log_message("Network visualization for {.val {pair}}")

  human_sample <- pair$human
  chimp_sample <- pair$chimp

  res_dir <- paste0(
    "results/species_networks/", human_sample, "_", chimp_sample, "/"
  )
  fig_dir <- check_dir(
    paste0("figures/species_networks/", human_sample, "_", chimp_sample)
  )

  node_annotation <- read.csv(
    file.path(res_dir, "network_evolution_core_table.csv"),
    stringsAsFactors = FALSE
  )

  trajectory_data <- read.csv(
    file.path(res_dir, "evolution_daccre_for_target_genes.csv"),
    stringsAsFactors = FALSE
  )

  csv_dir <- file.path(res_dir, "csv/")

  target_evolution_types <- c(
    "Regulatory rewiring",
    "Cis-regulatory activation",
    "Regulatory innovation",
    "Conserved"
  )

  for (selected_celltype in selected_celltypes) {
    if (!selected_celltype %in% node_annotation$CellType) {
      log_message(
        "Cell type {.val {selected_celltype}} not found in annotation, skipping...",
        message_type = "warning"
      )
      next
    }

    trajectory_ct <- trajectory_data[trajectory_data$CellType == selected_celltype, ]
    kept_targets <- trajectory_ct$Gene[
      trajectory_ct$Evolution_type %in% target_evolution_types
    ]

    edges_human_full <- read.csv(
      file.path(
        csv_dir,
        paste0("human_", human_sample, "_", selected_celltype, ".csv")
      ),
      stringsAsFactors = FALSE
    )
    edges_human_full <- edges_human_full[edges_human_full$target %in% kept_targets, ]
    edges_human_full <- edges_human_full[order(abs(edges_human_full$weight), decreasing = TRUE), ]
    targets_to_cover_h <- unique(edges_human_full$target)
    n_human <- if (length(targets_to_cover_h) == 0) {
      0L
    } else {
      covered <- character(0)
      for (i in seq_len(nrow(edges_human_full))) {
        covered <- union(covered, edges_human_full$target[i])
        if (all(targets_to_cover_h %in% covered)) break
      }
      i
    }
    log_message("Human targets to cover: {.val {length(targets_to_cover_h)}}")
    log_message("Human number of edges: {.val {n_human}}")
    edges_human <- head(edges_human_full, n_human)
    edges_human$edge_id <- paste(edges_human$regulator, edges_human$target, sep = "|")
    edges_human$source <- rep("Human", nrow(edges_human))

    chimp_file <- file.path(
      csv_dir,
      paste0("chimp_", chimp_sample, "_", selected_celltype, ".csv")
    )
    if (!file.exists(chimp_file)) {
      log_message(
        "Chimpanzee network file not found for {.val {selected_celltype}}, skipping...",
        message_type = "warning"
      )
      next
    }
    edges_chimp_full <- read.csv(chimp_file, stringsAsFactors = FALSE)
    edges_chimp_full <- edges_chimp_full[order(abs(edges_chimp_full$weight), decreasing = TRUE), ]
    targets_to_cover_c <- intersect(kept_targets, unique(edges_chimp_full$target))
    n_chimp <- if (length(targets_to_cover_c) == 0) {
      min(1000L, nrow(edges_chimp_full))
    } else {
      covered <- character(0)
      for (i in seq_len(nrow(edges_chimp_full))) {
        covered <- union(covered, edges_chimp_full$target[i])
        if (all(targets_to_cover_c %in% covered)) break
      }
      i
    }
    log_message("Chimpanzee targets to cover: {.val {length(targets_to_cover_c)}}")
    log_message("Chimpanzee number of edges: {.val {n_chimp}}")
    edges_chimp <- head(edges_chimp_full, n_chimp)
    edges_chimp$edge_id <- paste(
      edges_chimp$regulator, edges_chimp$target,
      sep = "|"
    )
    edges_chimp$source <- rep("Chimpanzee", nrow(edges_chimp))

    edges_merged <- merge(
      edges_human[, c("regulator", "target", "weight", "edge_id", "source")],
      edges_chimp[, c("regulator", "target", "weight", "edge_id", "source")],
      by = "edge_id",
      all = TRUE,
      suffixes = c("_human", "_chimp")
    )

    edges_merged$source <- ifelse(
      !is.na(edges_merged$source_human) & !is.na(edges_merged$source_chimp),
      "Not biased",
      ifelse(!is.na(edges_merged$source_human), "Human", "Chimpanzee")
    )

    edges_merged$weight <- ifelse(
      !is.na(edges_merged$weight_human),
      edges_merged$weight_human,
      edges_merged$weight_chimp
    )

    edges_merged$regulator <- ifelse(
      !is.na(edges_merged$regulator_human),
      edges_merged$regulator_human,
      edges_merged$regulator_chimp
    )
    edges_merged$target <- ifelse(
      !is.na(edges_merged$target_human),
      edges_merged$target_human,
      edges_merged$target_chimp
    )

    edges_merged <- edges_merged[, c("regulator", "target", "weight", "source")]

    g_merged <- graph_from_data_frame(
      edges_merged[, c("regulator", "target", "weight")],
      directed = TRUE,
      vertices = NULL
    )

    E(g_merged)$source <- edges_merged$source

    all_tfs <- unique(c(
      if (!is.null(edges_human)) unique(edges_human$regulator) else NULL,
      if (!is.null(edges_chimp)) unique(edges_chimp$regulator) else NULL
    ))

    V(g_merged)$type <- ifelse(V(g_merged)$name %in% all_tfs, "TF", "Target")

    node_anno_ct <- node_annotation[node_annotation$CellType == selected_celltype, ]
    V(g_merged)$network_status <- "Unknown"
    idx <- match(V(g_merged)$name, node_anno_ct$Gene)
    V(g_merged)$network_status[!is.na(idx)] <- node_anno_ct$Target_type[idx[!is.na(idx)]]
    V(g_merged)$network_status[V(g_merged)$type == "TF"] <- "TF"

    V(g_merged)$evolution_trajectory <- "Unknown"
    if (!is.null(trajectory_data)) {
      trajectory_ct <- trajectory_data[trajectory_data$CellType == selected_celltype, ]
      if (nrow(trajectory_ct) > 0) {
        idx_traj <- match(V(g_merged)$name, trajectory_ct$Gene)
        V(g_merged)$evolution_trajectory[!is.na(idx_traj)] <- trajectory_ct$Evolution_type[idx_traj[!is.na(idx_traj)]]
      }
    }

    V(g_merged)$evolution_trajectory[V(g_merged)$type == "TF"] <- "TF"

    V(g_merged)$evolution_trajectory[V(g_merged)$evolution_trajectory == "" | is.na(V(g_merged)$evolution_trajectory)] <- "Unknown"

    vert_remove <- V(g_merged)[
      V(g_merged)$type == "Target" &
        !(V(g_merged)$evolution_trajectory %in% target_evolution_types)
    ]
    g_merged <- igraph::delete_vertices(g_merged, vert_remove)

    ew <- abs(E(g_merged)$weight)
    V(g_merged)$strength <- igraph::strength(g_merged, weights = ew)
    innov_idx <- V(g_merged)$evolution_trajectory == "Regulatory innovation"
    if (sum(innov_idx) > 0) {
      strength_innov <- setNames(
        V(g_merged)$strength[innov_idx],
        V(g_merged)$name[innov_idx]
      )
      top10_innovation <- names(
        head(sort(strength_innov, decreasing = TRUE), top_n)
      )
    } else {
      top10_innovation <- character(0)
    }

    set.seed(2026)
    nv <- igraph::vcount(g_merged)
    layout_merged <- create_layout(
      g_merged,
      layout = "fr",
      weights = abs(E(g_merged)$weight),
      area = nv^2 * 4
    )

    p_merged <- ggraph(layout_merged) +
      geom_edge_link(
        aes(edge_width = abs(weight), edge_colour = source),
        alpha = 1,
        end_cap = circle(1, "mm")
      ) +
      scale_edge_width_continuous(range = c(0.1, 1), name = "Edge weight") +
      scale_edge_colour_manual(
        values = c(
          "Human" = human_color,
          "Chimpanzee" = chimp_color,
          "Not biased" = both_color
        ),
        name = "Edge type",
        guide = guide_legend(order = 2)
      ) +
      geom_node_point(
        aes(fill = evolution_trajectory, colour = "black", shape = type),
        size = 3,
        alpha = 1,
        stroke = 0.5
      ) +
      scale_fill_manual(
        values = trajectory_colors,
        name = "Evolution type",
        breaks = target_evolution_types,
        drop = FALSE,
        na.value = "#CCCCCC",
        guide = guide_legend(
          order = 1,
          override.aes = list(
            shape = 21,
            colour = "black",
            stroke = 0.5,
            size = 4
          )
        )
      ) +
      scale_colour_identity(guide = "none") +
      scale_shape_manual(
        values = c("TF" = 23, "Target" = 21),
        name = "Node type"
      ) +
      geom_node_text(
        data = function(x) x[x$name %in% top10_innovation, ],
        aes(label = name),
        repel = TRUE,
        size = 2.5,
        max.overlaps = 30,
        bg.colour = "white",
        bg.r = 0.15
      ) +
      theme_graph(base_family = "sans", background = "white") +
      theme(
        legend.position = "right",
        legend.box = "vertical",
        legend.key.size = unit(0.35, "cm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)
      )

    ggsave(
      file.path(
        fig_dir,
        paste0(
          "network_merged_", gsub(" ", "_", selected_celltype), ".pdf"
        )
      ),
      p_merged,
      width = 7,
      height = 5
    )

    log_message("Network plot saved for {.val {selected_celltype}}")
  }
}
