source("code/functions/aucell.R")
source("code/functions/hic.R")
source("code/functions/prepare_env.R")

edge_num <- 500
fig_dir <- check_dir("figures/aucell/")

intersection_results <- readRDS("results/hic/intersection_results.rds")

object <- readRDS(
  "../../data/BrainData/processed/GSE97942/GSE97942_cerebellum_processed.rds"
)

great_dataset <- read.csv("data/genome/science_HAR.csv")
great_dataset <- great_dataset |>
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
great_genes <- great_genes[great_genes %in% rownames(object)]

tfs_list <- extract_genes_list(
  intersection_results,
  type = "tfs",
  threshold = edge_num
)
tfs <- do.call(c, tfs_list) |> unique()

genes_list <- extract_genes_list(
  intersection_results,
  type = "all",
  threshold = edge_num
)
all_genes <- do.call(c, genes_list) |> unique()

gene_sets_list <- list(
  GENIE3 = genes_list$GENIE3,
  HARNexus = genes_list$HARNexus,
  LEAP = genes_list$LEAP,
  PPCOR = genes_list$PPCOR,
  GREAT = great_genes,
  TFs = tfs,
  All = all_genes
)

object <- aucell_analysis(
  seurat = object,
  highlight_celltypes = c("Astrocytes"),
  gene_sets_list = gene_sets_list,
  output_dir = paste0(fig_dir, edge_num, "_all"),
  threshold_offset = 0,
  use_adaptive_threshold = TRUE,
  adaptive_quantile_base = 0.9,
  adaptive_quantile_range = 0.1
)

genes_list <- extract_genes_list(
  intersection_results,
  type = "intersect",
  threshold = edge_num
)

gene_sets_list <- list(
  GENIE3 = genes_list$GENIE3,
  HARNexus = genes_list$HARNexus,
  LEAP = genes_list$LEAP,
  PPCOR = genes_list$PPCOR,
  GREAT = great_genes,
  TFs = tfs,
  All = all_genes
)

object <- aucell_analysis(
  seurat = object,
  highlight_celltypes = c("Astrocytes"),
  gene_sets_list = gene_sets_list,
  output_dir = paste0(fig_dir, edge_num, "_intersect"),
  threshold_offset = 0,
  use_adaptive_threshold = TRUE,
  adaptive_quantile_base = 0.9,
  adaptive_quantile_range = 0.1
)

genes_list <- extract_genes_list(
  intersection_results,
  type = "setdiff",
  threshold = edge_num
)

gene_sets_list <- list(
  GENIE3 = genes_list$GENIE3,
  HARNexus = genes_list$HARNexus,
  LEAP = genes_list$LEAP,
  PPCOR = genes_list$PPCOR,
  GREAT = great_genes,
  TFs = tfs,
  All = all_genes
)

object <- aucell_analysis(
  seurat = object,
  highlight_celltypes = c("Astrocytes"),
  gene_sets_list = gene_sets_list,
  output_dir = paste0(fig_dir, edge_num, "_setdiff"),
  threshold_offset = 0,
  use_adaptive_threshold = TRUE,
  adaptive_quantile_base = 0.9,
  adaptive_quantile_range = 0.1
)
