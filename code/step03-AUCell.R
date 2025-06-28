source("AUCell.R")
source("functions_hic.R")

data_dir <- "results/hi-c/data"
figures_paper <- "results/plots_paper"
figures_zhang <- "results/plots"

load(paste0(data_dir, "/input_data.Rdata"))
load(paste0(data_dir, "/intersection_results.Rdata"))

tfs <- read.table("data/regulators.txt", header = FALSE)[, 1]

if (!file.exists("data/seurat_sub.Rdata")) {
  load("data/CBC_S8_seurat.Rdata")

  cell_counts <- table(seurat$cell_name)
  major_cell_types <- names(cell_counts[cell_counts > 20])

  seurat <- subset(seurat, cell_name %in% major_cell_types)
  print(table(seurat$cell_name))
  seurat <- seurat_functions(seurat)

  seurat_sub <- subset(seurat, cell_name == "Astro")
  save(seurat, seurat_sub, file = "data/seurat_sub.Rdata")
} else {
  load("data/seurat_sub.Rdata")
}
seurat <- Seurat::RunUMAP(
  seurat,
  dims = 1:10,
  reduction = "pca",
  n.neighbors = 30,
  min.dist = 0.3
)
seurat <- Seurat::RunTSNE(
  seurat,
  tsne.method = "Rtsne",
  dim.embed = 2,
  reduction.key = "tSNE_"
)
expr_matrix <- as.matrix(
  GetAssayData(seurat, layer = "data")
)

great_dataset <- readr::read_csv("data/science_HAR.csv")
great_dataset <- great_dataset %>%
  as_tibble() %>%
  separate_rows(associated_genes, sep = ",")
great_dataset$TSS <- gsub(
  ".*\\(|\\).*", "",
  great_dataset$associated_genes
)
great_dataset$associated_genes <- gsub(
  "\\(.*?\\)", "",
  great_dataset$associated_genes
)
great_dataset$associated_genes <- gsub(
  " ", "",
  great_dataset$associated_genes
)
great_genes <- unique(great_dataset$associated_genes)
great_genes <- great_genes[great_genes %in% rownames(expr_matrix)]
length(great_genes)

genes_list <- extract_genes_list(
  intersection_results,
  type = "all",
  threshold = 1000
)

gene_sets_list <- list(
  GENIE3 = genes_list$GENIE3,
  LEAP = genes_list$LEAP,
  PPCOR = genes_list$PPCOR,
  HARNexus = genes_list$HARNexus,
  GREAT = great_genes
)

results <- aucell_analysis(
  seurat = seurat,
  gene_sets_list = gene_sets_list,
  random_sets = TRUE,
  random_sizes = c(100, 200, 400),
  output_dir = paste0(figures_paper, "/AUCell_1000_all"),
  seed = 2024,
  threshold_offset = 0.02,
  point_size = 0.5
)


genes_list <- extract_genes_list(
  intersection_results,
  type = "intersect",
  threshold = 1000
)

gene_sets_list <- list(
  GENIE3 = genes_list$GENIE3,
  LEAP = genes_list$LEAP,
  PPCOR = genes_list$PPCOR,
  HARNexus = genes_list$HARNexus,
  GREAT = great_genes
)

results <- aucell_analysis(
  seurat = seurat,
  gene_sets_list = gene_sets_list,
  random_sets = TRUE,
  random_sizes = c(100, 200, 400),
  output_dir = paste0(figures_paper, "/AUCell_1000_intersect"),
  threshold_offset = 0.01,
  point_size = 0.5
)

genes_list <- extract_genes_list(
  intersection_results,
  type = "setdiff",
  threshold = 1000
)

gene_sets_list <- list(
  GENIE3 = genes_list$GENIE3,
  LEAP = genes_list$LEAP,
  PPCOR = genes_list$PPCOR,
  HARNexus = genes_list$HARNexus,
  GREAT = great_genes
)

results <- aucell_analysis(
  seurat = seurat,
  gene_sets_list = gene_sets_list,
  random_sets = TRUE,
  random_sizes = c(100, 200, 400),
  output_dir = paste0(figures_paper, "/AUCell_1000_setdiff"),
  seed = 2024,
  threshold_offset = 0.02,
  point_size = 0.5
)

tfs_infercsn <- unique(network_table_infercsn$regulator)
targets_infercsn <- unique(network_table_infercsn$target)

network_table_genie3_infercsn <- network_table_genie3[network_table_genie3$regulator %in% tfs_infercsn, ]
network_table_genie3_infercsn <- network_table_genie3_infercsn[network_table_genie3_infercsn$target %in% targets_infercsn, ]
network_table_leap_infercsn <- network_table_leap[network_table_leap$regulator %in% tfs_infercsn, ]
network_table_leap_infercsn <- network_table_leap_infercsn[network_table_leap_infercsn$target %in% targets_infercsn, ]
network_table_ppcor_infercsn <- network_table_ppcor[network_table_ppcor$regulator %in% tfs_infercsn, ]
network_table_ppcor_infercsn <- network_table_ppcor_infercsn[network_table_ppcor_infercsn$target %in% targets_infercsn, ]

genes_list_infercsn <- list(
  GENIE3 = as.character(unique(network_table_genie3_infercsn$target)),
  LEAP = as.character(unique(network_table_leap_infercsn$target)),
  PPCOR = as.character(unique(network_table_ppcor_infercsn$target)),
  HARNexus = as.character(targets_infercsn),
  GREAT = as.character(great_genes)
)

results <- aucell_analysis(
  seurat = seurat,
  gene_sets_list = genes_list_infercsn,
  random_sets = TRUE,
  random_sizes = c(100, 200, 400),
  output_dir = paste0(figures_paper, "/AUCell_1000_infercsn"),
  seed = 2024,
  threshold_offset = 0.02,
  point_size = 0.5
)
