get_celltype_specific_genes <- function(
    seurat,
    celltype_col = "CellType",
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
    stop("Input 'seurat' must be a Seurat object.")
  }

  if (!celltype_col %in% colnames(seurat@meta.data)) {
    stop(
      "Cell type column '", celltype_col,
      "' not found in Seurat object metadata."
    )
  }

  if (verbose) {
    log_message("Extracting expression matrix from assay: {.val {assay}}, layer: {.val {layer}}")
  }
  expr_mat <- Seurat::GetAssayData(seurat, assay = assay, layer = layer)

  if (nrow(expr_mat) == 0 || ncol(expr_mat) == 0) {
    stop("Expression matrix is empty.")
  }

  cell_type_vec <- seurat@meta.data[[celltype_col]]
  if (is.null(cell_type_vec)) {
    stop("Failed to extract cell type information from metadata.")
  }

  cell_types <- sort(unique(cell_type_vec))

  if (length(cell_types) == 0) {
    stop("No cell types found in the specified column.")
  }

  if (verbose) {
    log_message("Found {.val {length(cell_types)}} cell types")
  }

  gene_sets_list <- list()

  for (ct in cell_types) {
    if (verbose) {
      log_message("Processing cell type: {.val {ct}}")
    }

    cells_ct_idx <- which(cell_type_vec == ct)
    n_ct <- length(cells_ct_idx)

    if (n_ct == 0) {
      warning("No cells found for cell type: ", ct, ", skipping.")
      next
    }

    if (n_ct < min_cells) {
      min_prop_expressed <- min(0.05, min_prop_expressed_base)
    } else {
      min_prop_expressed <- min_prop_expressed_base
    }

    expr_ct <- expr_mat[, cells_ct_idx, drop = FALSE]

    prop_expr <- Matrix::rowMeans(expr_ct > 0)
    genes_prop_keep_idx <- which(prop_expr >= min_prop_expressed)

    if (verbose) {
      log_message("  After prop filter: {.val {length(genes_prop_keep_idx)}} genes")
    }

    if (length(genes_prop_keep_idx) == 0) {
      log_message(
        "  No genes passed prop filter for {.val {ct}}, skip.",
        message_type = "warning"
      )
      next
    }

    mean_expr_ct_all <- Matrix::rowMeans(expr_ct)
    mean_expr_ct_sel <- mean_expr_ct_all[genes_prop_keep_idx]

    genes_mean_keep_idx <- genes_prop_keep_idx[mean_expr_ct_sel >= min_mean_expr]

    if (verbose) {
      log_message("  After mean filter: {.val {length(genes_mean_keep_idx)}} genes")
    }

    if (length(genes_mean_keep_idx) == 0) {
      log_message(
        "  No genes passed mean filter for {.val {ct}}, skip.",
        message_type = "warning"
      )
      next
    }

    expr_ct_dense <- as.matrix(expr_ct[genes_mean_keep_idx, , drop = FALSE])
    sd_expr_ct <- matrixStats::rowSds(expr_ct_dense)

    genes_sd_keep_idx <- genes_mean_keep_idx[sd_expr_ct >= min_sd_expr]

    if (verbose) {
      log_message("  After SD filter: {.val {length(genes_sd_keep_idx)}} genes")
    }

    if (length(genes_sd_keep_idx) == 0) {
      log_message(
        "  No genes passed SD filter for {.val {ct}}, skip.",
        message_type = "warning"
      )
      next
    }

    # Filter 4: Specificity Index (SI)
    # SI_g,T = mu_g,T / (max_{others} mu_g,others + epsilon)
    expr_ct_all <- expr_mat[, cells_ct_idx, drop = FALSE]
    mu_ct_all <- Matrix::rowMeans(expr_ct_all)

    cells_other_idx <- which(cell_type_vec != ct)
    if (length(cells_other_idx) > 0) {
      expr_other <- expr_mat[, cells_other_idx, drop = FALSE]
      mu_other_all <- Matrix::rowMeans(expr_other)
    } else {
      # If no other cells, set mu_other_all to zero
      mu_other_all <- rep(0, nrow(expr_mat))
      names(mu_other_all) <- rownames(expr_mat)
    }

    SI_all <- mu_ct_all / (mu_other_all + epsilon)

    SI_sel <- SI_all[genes_sd_keep_idx]
    genes_SI_keep_idx <- genes_sd_keep_idx[SI_sel >= min_SI]

    if (verbose) {
      log_message("  After SI filter: {.val {length(genes_SI_keep_idx)}} genes")
    }

    if (length(genes_SI_keep_idx) == 0) {
      warning("  No genes passed SI filter for ", ct, ", skip.")
      next
    }

    # Extract final gene names
    final_genes_ct <- rownames(expr_mat)[genes_SI_keep_idx]
    gene_sets_list[[ct]] <- final_genes_ct

    if (verbose) {
      log_message("  Final genes for {.val {ct}}: {.val {length(final_genes_ct)}}")
    }
  }

  if (length(gene_sets_list) == 0) {
    warning("No cell type-specific genes found for any cell type.")
  }

  return(gene_sets_list)
}
