source("code/functions/prepare_env.R")

sample_pairs <- list(
  list(human = "h4", chimp = "c4"),
  list(human = "h3", chimp = "c2"),
  list(human = "h1", chimp = "c1")
)

for (pair in sample_pairs) {
  human_sample <- pair$human
  chimp_sample <- pair$chimp

  log_message(
    "DAcCRE-based analysis for {.val {pair}}"
  )

  res_dir <- paste0(
    "results/species_networks/", human_sample, "_", chimp_sample, "/"
  )
  atac_dir <- paste0(
    "results/species_networks/", human_sample, "_", chimp_sample, "/atac"
  )

  genes_human <- read.csv(
    file.path(res_dir, "network_genes_by_celltype_human.csv"),
    stringsAsFactors = FALSE
  )
  genes_chimp <- read.csv(
    file.path(res_dir, "network_genes_by_celltype_chimp.csv"),
    stringsAsFactors = FALSE
  )

  targets_human <- genes_human[!genes_human$is_TF, ]
  targets_chimp <- genes_chimp[!genes_chimp$is_TF, ]

  log_message("Human targets: {.val {length(unique(targets_human$gene))}}")
  log_message("Chimp targets: {.val {length(unique(targets_chimp$gene))}}")

  all_celltypes <- union(
    unique(targets_human$CellType),
    unique(targets_chimp$CellType)
  )
  log_message("Cell types (Human or Chimp): {.val {all_celltypes}}")

  results <- list()
  for (ct in all_celltypes) {
    human_genes_ct <- unique(targets_human$gene[targets_human$CellType == ct])
    chimp_genes_ct <- unique(targets_chimp$gene[targets_chimp$CellType == ct])
    all_genes_ct <- unique(c(human_genes_ct, chimp_genes_ct))

    for (g in all_genes_ct) {
      in_human <- g %in% human_genes_ct
      in_chimp <- g %in% chimp_genes_ct
      network_status <- ifelse(
        in_human & !in_chimp, "Human-only",
        ifelse(!in_human & in_chimp, "Chimp-only", "Not biased")
      )
      results[[paste(ct, g, sep = "_")]] <- data.frame(
        Gene = g,
        CellType = ct,
        Target_type = network_status,
        stringsAsFactors = FALSE
      )
    }
  }

  core_df <- do.call(rbind, results)
  write.csv(
    core_df,
    file.path(res_dir, "network_evolution_core_table.csv"),
    row.names = FALSE
  )

  p2g <- read.csv(
    file.path(atac_dir, "peak2gene_500kb.csv"),
    stringsAsFactors = FALSE
  )
  gene_peaks <- split(p2g$peak_id, p2g$gene)

  daccre_files <- list.files(
    atac_dir,
    pattern = "^daccre_.*\\.csv$", full.names = TRUE
  )

  peak_status_by_ct <- list()
  for (f in daccre_files) {
    ct <- gsub("^daccre_", "", basename(f))
    ct <- gsub("\\.csv$", "", ct)
    ct <- gsub("_", " ", ct)
    if (ct == "Oligodendrocyte") {
      next
    }
    d <- read.csv(f, stringsAsFactors = FALSE)
    peak_status_by_ct[[ct]] <- d[, c("peak_id", "peak_status")]
  }

  core_df$Peak_type <- NA_character_
  core_df$n_peaks <- 0L
  core_df$n_Human_gained <- 0L
  core_df$n_Chimp_gained <- 0L
  core_df$n_both_peaks <- 0L

  for (i in seq_len(nrow(core_df))) {
    g <- core_df$Gene[i]
    ct <- core_df$CellType[i]
    linked_peaks <- gene_peaks[[g]]

    if (is.null(linked_peaks) || length(linked_peaks) == 0) {
      core_df$Peak_type[i] <- "Chromatin inaccessible"
      next
    }

    peak_status_df <- peak_status_by_ct[[ct]]
    if (is.null(peak_status_df)) {
      core_df$Peak_type[i] <- "No_data"
      next
    }

    matched_peaks <- peak_status_df[peak_status_df$peak_id %in% linked_peaks, ]
    if (nrow(matched_peaks) == 0) {
      core_df$Peak_type[i] <- "Chromatin inaccessible"
      next
    }

    core_df$n_peaks[i] <- nrow(matched_peaks)
    core_df$n_Human_gained[i] <- sum(
      matched_peaks$peak_status == "Human-biased"
    )
    core_df$n_Chimp_gained[i] <- sum(
      matched_peaks$peak_status == "Chimp-biased"
    )
    is_shared_peak <- matched_peaks$peak_status %in% "Not-biased"
    core_df$n_both_peaks[i] <- sum(is_shared_peak)
    has_shared_peak <- any(is_shared_peak)

    has_hg <- any(matched_peaks$peak_status == "Human-biased")
    has_cg <- any(matched_peaks$peak_status == "Chimp-biased")

    # All branches are mutually exclusive:
    # - Chromatin accessible (human): only Human-biased peaks
    # - Chromatin accessible (chimp): only Chimp-biased peaks
    # - Chromatin inaccessible: no peaks
    if (has_hg && !has_cg) {
      core_df$Peak_type[i] <- "Chromatin accessible"
    } else if (has_cg && !has_hg) {
      core_df$Peak_type[i] <- "Chromatin inaccessible"
    } else if (has_shared_peak && !has_hg && !has_cg) {
      core_df$Peak_type[i] <- "Chromatin accessible"
    } else {
      core_df$Peak_type[i] <- "Chromatin inaccessible"
    }
  }

  core_df$Evolution_type <- NA_character_

  for (i in seq_len(nrow(core_df))) {
    net_status <- core_df$Target_type[i]
    peak_status <- core_df$Peak_type[i]

    if (is.na(net_status) || is.na(peak_status)) {
      core_df$Evolution_type[i] <- "Other"
      next
    }

    # Create evolution trajectory classification
    # - Regulatory Innovation: new element drives new network
    # - Cis-regulatory Activation: old connection, new activity
    # - Regulatory Rewiring: old activity, new connection
    # - Conserved: new network and old activity
    if (net_status == "Human-only" && peak_status == "Chromatin accessible") {
      core_df$Evolution_type[i] <- "Regulatory innovation"
    } else if (net_status == "Not biased" && peak_status == "Chromatin accessible") {
      core_df$Evolution_type[i] <- "Cis-regulatory activation"
    } else if (net_status == "Human-only" && peak_status == "Chromatin inaccessible") {
      core_df$Evolution_type[i] <- "Regulatory rewiring"
    } else {
      core_df$Evolution_type[i] <- "Conserved"
    }
  }

  log_message("Evolution types assigned")

  write.csv(
    core_df,
    file.path(res_dir, "evolution_daccre_for_target_genes.csv"),
    row.names = FALSE
  )

  trajectory_summary <- as.data.frame(
    table(
      Target_type = core_df$Target_type,
      Peak_type = core_df$Peak_type,
      Evolution_type = core_df$Evolution_type
    )
  )
  colnames(trajectory_summary)[4] <- "Count"

  trajectory_summary <- trajectory_summary[trajectory_summary$Count > 0, ]
  trajectory_summary <- trajectory_summary[order(trajectory_summary$Count, decreasing = TRUE), ]

  write.csv(
    trajectory_summary,
    file.path(res_dir, "evolution_trajectory_summary.csv"),
    row.names = FALSE
  )

  trajectory_by_ct_list <- list()
  celltypes <- unique(core_df$CellType)

  for (ct in celltypes) {
    ct_data <- core_df[core_df$CellType == ct, ]
    traj_table <- table(ct_data$Evolution_type)

    traj_df <- data.frame(
      CellType = ct,
      Evolution_type = names(traj_table),
      Count = as.numeric(traj_table),
      Proportion = round(as.numeric(traj_table) / nrow(ct_data), 3),
      stringsAsFactors = FALSE
    )
    traj_df <- traj_df[order(traj_df$Count, decreasing = TRUE), ]
    trajectory_by_ct_list[[length(trajectory_by_ct_list) + 1]] <- traj_df
  }

  trajectory_by_ct_flat <- do.call(rbind, trajectory_by_ct_list)

  write.csv(
    trajectory_by_ct_flat,
    file.path(res_dir, "evolution_types_by_celltype.csv"),
    row.names = FALSE
  )

  key_trajectories <- c(
    "Regulatory innovation",
    "Cis-regulatory activation",
    "Regulatory rewiring",
    "Conserved"
  )

  key_summary <- data.frame(
    Evolution_type = key_trajectories,
    Count = sapply(key_trajectories, function(t) {
      sum(core_df$Evolution_type == t, na.rm = TRUE)
    }),
    stringsAsFactors = FALSE
  )
  key_summary$Proportion <- round(key_summary$Count / nrow(core_df), 3)

  write.csv(
    key_summary,
    file.path(res_dir, "evolution_types_summary.csv"),
    row.names = FALSE
  )
}
