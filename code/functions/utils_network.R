get_celltype_specific_genes <- function(
    seurat,
    group.by = "CellType",
    assay = "RNA",
    layer = "data",
    min_prop_expressed_base = 0.10,
    min_mean_expr = 0.10,
    min_sd_expr = 0.20,
    min_SI = 1.20,
    epsilon = 1e-6,
    min_cells = 100,
    verbose = TRUE) {
  if (!inherits(seurat, "Seurat")) {
    log_message(
      "Input 'seurat' must be a Seurat object",
      message_type = "error"
    )
  }

  if (!group.by %in% colnames(seurat@meta.data)) {
    log_message(
      "{.val {group.by}} not found in {.cls Seurat}",
      message_type = "error"
    )
  }

  log_message(
    "Extracting expression matrix from assay: {.val {assay}}, layer: {.val {layer}}",
    verbose = verbose
  )
  expr_data <- Seurat::GetAssayData(seurat, assay = assay, layer = layer)
  expr_data <- expr_data[rowSums(expr_data) > 0, ]
  if (nrow(expr_data) == 0 || ncol(expr_data) == 0) {
    log_message(
      "Expression matrix is empty",
      message_type = "error"
    )
  }

  cell_type_vec <- seurat@meta.data[[group.by]]
  if (is.null(cell_type_vec)) {
    log_message(
      "Failed to extract cell type information from metadata",
      message_type = "error"
    )
  }

  cell_types <- sort(unique(cell_type_vec))

  if (length(cell_types) == 0) {
    log_message(
      "No cell types found in the specified column",
      message_type = "error"
    )
  }

  log_message(
    "Found {.val {length(cell_types)}} cell types",
    verbose = verbose
  )

  gene_sets_list <- list()

  for (ct in cell_types) {
    log_message("Processing cell type: {.val {ct}}", verbose = verbose)

    cells_ct_idx <- which(cell_type_vec == ct)
    n_ct <- length(cells_ct_idx)

    if (n_ct == 0) {
      log_message(
        "No cells found for cell type: {.val {ct}}, skipping",
        message_type = "warning",
        verbose = verbose
      )
      next
    }

    if (n_ct < min_cells) {
      min_prop_expressed <- min(0.05, min_prop_expressed_base)
    } else {
      min_prop_expressed <- min_prop_expressed_base
    }

    expr_ct <- expr_data[, cells_ct_idx, drop = FALSE]

    prop_expr <- Matrix::rowMeans(expr_ct > 0)
    genes_prop_keep_idx <- which(prop_expr >= min_prop_expressed)

    if (length(genes_prop_keep_idx) == 0) {
      log_message(
        "  No genes passed prop filter for {.val {ct}}, skip",
        message_type = "warning",
        verbose = verbose
      )
      next
    }

    mean_expr_ct_all <- Matrix::rowMeans(expr_ct)
    mean_expr_ct_sel <- mean_expr_ct_all[genes_prop_keep_idx]

    genes_mean_keep_idx <- genes_prop_keep_idx[mean_expr_ct_sel >= min_mean_expr]

    if (length(genes_mean_keep_idx) == 0) {
      log_message(
        "  No genes passed mean filter for {.val {ct}}, skip",
        message_type = "warning",
        verbose = verbose
      )
      next
    }

    expr_ct_dense <- as.matrix(expr_ct[genes_mean_keep_idx, , drop = FALSE])
    sd_expr_ct <- matrixStats::rowSds(expr_ct_dense)

    genes_sd_keep_idx <- genes_mean_keep_idx[sd_expr_ct >= min_sd_expr]

    if (length(genes_sd_keep_idx) == 0) {
      log_message(
        "  No genes passed SD filter for {.val {ct}}, skip.",
        message_type = "warning",
        verbose = verbose
      )
      next
    }

    # Filter 4: Specificity Index (SI)
    # SI_g,T = mu_g,T / (max_{others} mu_g,others + epsilon)
    expr_ct_all <- expr_data[, cells_ct_idx, drop = FALSE]
    mu_ct_all <- Matrix::rowMeans(expr_ct_all)

    cells_other_idx <- which(cell_type_vec != ct)
    if (length(cells_other_idx) > 0) {
      expr_other <- expr_data[, cells_other_idx, drop = FALSE]
      mu_other_all <- Matrix::rowMeans(expr_other)
    } else {
      mu_other_all <- rep(0, nrow(expr_data))
      names(mu_other_all) <- rownames(expr_data)
    }

    SI_all <- mu_ct_all / (mu_other_all + epsilon)

    SI_sel <- SI_all[genes_sd_keep_idx]
    genes_SI_keep_idx <- genes_sd_keep_idx[SI_sel >= min_SI]

    if (length(genes_SI_keep_idx) == 0) {
      log_message(
        "  No genes passed SI filter for {.val {ct}}, skip",
        message_type = "warning",
        verbose = verbose
      )
      next
    }

    final_genes_ct <- rownames(expr_data)[genes_SI_keep_idx]
    gene_sets_list[[ct]] <- final_genes_ct

    log_message(
      "  Final genes for {.val {ct}}: {.val {length(final_genes_ct)}}",
      verbose = verbose
    )
  }

  if (length(gene_sets_list) == 0) {
    log_message(
      "No cell type-specific genes found for any cell type",
      message_type = "warning",
      verbose = verbose
    )
  }

  return(gene_sets_list)
}


run_genie3 <- function(
    matrix,
    regulators = NULL,
    targets = NULL,
    cores = 1) {
  network_table <- GENIE3::getLinkList(
    GENIE3::GENIE3(
      exprMatrix = as.matrix(matrix),
      regulators = regulators,
      targets = targets,
      verbose = TRUE,
      nCores = cores
    )
  )
  colnames(network_table) <- c("regulator", "target", "weight")
  return(network_table)
}

run_ppcor <- function(
    matrix,
    regulators = NULL,
    targets = NULL) {
  gene_names <- rownames(matrix)
  rownames(matrix) <- gene_names
  ppcor_results <- ppcor::pcor(
    x = t(as.matrix(matrix)),
    method = "spearman"
  )
  network_table <- data.frame(
    regulator = gene_names[c(row(ppcor_results$estimate))],
    target = gene_names[c(col(ppcor_results$estimate))],
    weight = c(ppcor_results$estimate),
    p_value = c(ppcor_results$p.value)
  )

  if (!is.null(regulators)) {
    network_table <- network_table[network_table$regulator %in% regulators, ]
  }
  if (!is.null(targets)) {
    network_table <- network_table[network_table$target %in% targets, ]
  }
  network_table <- network_table[network_table$weight != 0, ]

  network_table$weight <- as.numeric(network_table$weight)

  network_table <- network_table[order(
    abs(network_table$weight),
    decreasing = TRUE
  ), ][, -4]

  colnames(network_table) <- c("regulator", "target", "weight")
  return(network_table)
}

run_leap <- function(
    matrix,
    regulators = NULL,
    targets = NULL) {
  gene_names <- rownames(matrix)
  rownames(matrix) <- c()
  leap_results <- LEAP::MAC_counter(
    data = matrix,
    max_lag_prop = 1 / 3,
    MAC_cutoff = 0.2,
    file_name = FALSE,
    lag_matrix = TRUE,
    symmetric = FALSE
  )

  network_table <- data.frame(
    regulator = gene_names[leap_results[, "Row gene index"]],
    target = gene_names[leap_results[, "Column gene index"]],
    weight = leap_results[, "Correlation"]
  )
  if (!is.null(regulators)) {
    network_table <- network_table[network_table$regulator %in% regulators, ]
  }
  if (!is.null(targets)) {
    network_table <- network_table[network_table$target %in% targets, ]
  }
  network_table <- network_table[network_table$weight != 0, ]

  network_table$weight <- as.numeric(network_table$weight)
  network_table <- network_table[order(
    abs(network_table$weight),
    decreasing = TRUE
  ), ]

  colnames(network_table) <- c("regulator", "target", "weight")
  return(network_table)
}
