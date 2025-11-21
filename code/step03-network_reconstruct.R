source("functions.R")

tfs <- read.table("data/regulators.txt", header = FALSE)[, 1]
load("data/CBC_S8_seurat.Rdata")
table(seurat$cluster)

cell_counts <- table(seurat$cluster)
major_cell_types <- names(cell_counts[cell_counts > 50])

seurat <- subset(seurat, cluster %in% major_cell_types)
print(table(seurat$cluster))
seurat <- seurat_functions(seurat)
expr_matrix <- as.matrix(
  GetAssayData(seurat, layer = "data")
)
var_genes <- VariableFeatures(seurat)
seurat_sub <- subset(seurat, cluster == "Astro")

expr_matrix_astro <- as.matrix(
  GetAssayData(seurat_sub, layer = "data")
)
expr_matrix_astro <- expr_matrix_astro[var_genes, ]
expr_matrix_astro <- expr_matrix_astro[rowSums(expr_matrix_astro) > 0, ]
regulators <- tfs[tfs %in% rownames(expr_matrix_astro)]

if (!file.exists("data/Astro_network_table_list.Rdata")) {
  network_table_genie3 <- run_GENIE3(
    expr_matrix_astro,
    regulators = regulators
  )
  network_table_ppcor <- run_ppcor(
    expr_matrix_astro,
    regulators = regulators
  )
  network_table_leap <- run_LEAP(
    expr_matrix_astro,
    regulators = regulators
  )

  load("data/CBC_S8_grn_inferCSN.Rdata")
  network_table_infercsn <- grn_list$Astro
  save(
    network_table_genie3,
    network_table_ppcor,
    network_table_leap,
    network_table_infercsn,
    file = "data/Astro_network_table_list.Rdata"
  )
}

load("data/Astro_network_table_list.Rdata")
