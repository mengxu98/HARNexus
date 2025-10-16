rm(list = ls())
gc()

source("code/functions/packages.R")

# paper: https://doi.org/10.1016/j.xgen.2024.100703
# code: https://github.com/KIZ-SubLab/ACC-sn-Multiomes
# data: https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA015229

data_dir <- "../../data/BrainData/"
res_dir <- "../../data/BrainData/results/2024CellGenomics/"

log_message("Start loading data...")

data_path <- read.csv(
  file.path(data_dir, "2024CellGenomics/PRJCA015229/data.detail.csv")
)

data_list <- list()
data_path_multiomes <- subset(data_path, data.type == "multiomes")
data_path_multiomes <- subset(data_path_multiomes, species == "Human")
data_path_rna <- subset(data_path, data.type == "snRNA")
data_path_rna <- subset(data_path_rna, species == "Human")

for (i in seq_along(data_path_multiomes$objectID)) {
  data_list[[i]] <- Read10X(
    file.path(
      data_dir,
      gsub("data/", "", data_path_multiomes$pathway[i]),
      "filtered_feature_bc_matrix"
    ),
    gene.column = 2
  )
  data_list[[i]] <- data_list[[i]][[1]]
}
for (i in seq_along(data_path_rna$objectID)) {
  data_list[[i + length(data_path_multiomes$objectID)]] <- Read10X(
    file.path(
      data_dir,
      gsub("data/", "", data_path_rna$pathway[i]),
      "filtered_feature_bc_matrix"
    ),
    gene.column = 2
  )
}

human_id <- purrr::map(
  data_list, function(x) unique(rownames(x))
) |>
  purrr::list_c() |>
  unique()

for (i in seq_along(data_list)) {
  data_list[[i]] <- data_list[[i]][human_id, ]
}
names(data_list) <- c(data_path_multiomes$objectID, data_path_rna$objectID)

object_list <- list()
for (i in seq_along(data_list)) {
  object_list[[i]] <- CreateSeuratObject(
    counts = data_list[[i]],
    project = names(data_list)[i]
  )
  object_list[[i]][["homo"]] <- CreateAssayObject(
    counts = data_list[[i]]
  )
}
names(object_list) <- names(data_list)
rm(data_list)
gc()

saveRDS(
  object_list,
  file.path(
    res_dir, "object_list(6samples)_snRNA-seq.rds"
  )
)

object_list <- readRDS(
  file.path(
    res_dir, "object_list(6samples)_snRNA-seq.rds"
  )
)
for (i in seq_along(object_list)) {
  object_list[[i]] <- scop::RunDoubletCalling(
    object_list[[i]],
    db_method = "Scrublet"
  )
  object_list[[i]]@meta.data[["DoubletScores"]] <- object_list[[i]]@meta.data[["db.Scrublet_score"]]
  object_list[[i]]@meta.data[["PredictedDoublets"]] <- object_list[[i]]@meta.data[["db.Scrublet_class"]]
  object_list[[i]]@meta.data[["DoubletScores"]] <- unlist(object_list[[i]]@meta.data[["DoubletScores"]])
  object_list[[i]]@meta.data[["PredictedDoublets"]] <- unlist(object_list[[i]]@meta.data[["PredictedDoublets"]])
}

for (i in seq_along(object_list)) {
  object_list[[i]][["percent_mito"]] <- PercentageFeatureSet(
    object_list[[i]],
    pattern = "^MT-"
  )
}
object_list <- lapply(
  X = object_list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(
      x,
      selection.method = "vst", nfeatures = 2000, verbose = FALSE
    )
    x <- ScaleData(x, verbose = FALSE)
    x <- RunPCA(x, npcs = 30, verbose = FALSE)
    x <- RunUMAP(x, reduction = "pca", dims = 1:30, verbose = FALSE)
    x <- FindNeighbors(x, reduction = "pca", dims = 1:30, verbose = FALSE)
    x <- FindClusters(x, resolution = 0.5, verbose = FALSE)
  }
)

saveRDS(
  object_list,
  file.path(
    res_dir, "object_list(6samples)_snRNA-seq.rds"
  )
)

object_list <- readRDS(
  file.path(
    res_dir, "object_list(6samples)_snRNA-seq.rds"
  )
)

# object_filter_list <- list()
# for (i in seq_along(object_list)) {
#   object_filter_list[[i]] <- subset(
#     object_list[[i]],
#     subset =
#       nFeature_RNA > 300 &
#         nFeature_RNA < 5000 &
#         percent_mito < 10 &
#         DoubletScores < 1 &
#         PredictedDoublets != "doublets"
#   )
# }
# rm(object_list)
# gc()

# log_message("Finding common genes across all samples")
# log_message("Genes per sample before intersection:")
# for (i in seq_along(object_filter_list)) {
#   log_message("Sample", i, ":", nrow(object_filter_list[[i]]), "genes")
# }
# common_genes <- rownames(object_filter_list[[1]])
# for (i in 2:length(object_filter_list)) {
#   common_genes <- intersect(common_genes, rownames(object_filter_list[[i]]))
# }
# log_message("Number of common genes:", length(common_genes))

# for (i in seq_along(object_filter_list)) {
#   object_filter_list[[i]] <- subset(
#     object_filter_list[[i]],
#     features = common_genes
#   )
# }

# for (i in seq_along(object_filter_list)) {
#   DefaultAssay(object_filter_list[[i]]) <- "homo"
# }
# object_filter_list <- lapply(
#   X = object_filter_list,
#   FUN = SCTransform,
#   method = "glmGamPoi"
# )
# features <- SelectIntegrationFeatures(
#   object.list = object_filter_list,
#   nfeatures = 2000
# )
# object_filter_list <- PrepSCTIntegration(
#   object.list = object_filter_list,
#   anchor.features = features
# )
# object_filter_list <- lapply(
#   X = object_filter_list,
#   FUN = RunPCA,
#   features = features
# )
# anchors <- FindIntegrationAnchors(
#   object.list = object_filter_list,
#   normalization.method = "SCT",
#   anchor.features = features,
#   reduction = "cca",
#   dims = 1:30
# )

# log_message("integrate data")
# object <- IntegrateData(
#   anchorset = anchors,
#   normalization.method = "SCT",
#   dims = 1:30
# )

# DefaultAssay(object) <- "RNA"


object <- merge(object_list[[1]], object_list[2:6])
object <- JoinLayers(object)

object <- subset(
  object,
  PredictedDoublets == "singlet"
)

srt_cells <- colnames(object)
cell_ids <- gsub("_[0-9]+$", "", srt_cells)
orig_ident <- object$orig.ident
new_colnames <- paste0(cell_ids, "_", orig_ident)
colnames(object) <- new_colnames

# https://www.cell.com/cms/10.1016/j.xgen.2024.100703/attachment/29b2bce9-946b-43b2-8f12-346dbaea1b6f/mmc6.xlsx
metadata_snmultiome <- read.csv(
  file.path(data_dir, "2024CellGenomics/metadata_snMultiome.csv")
)
metadata_snrnaseq <- read.csv(
  file.path(data_dir, "2024CellGenomics/metadata_snRNA-seq.csv")
)
common_colnames <- intersect(
  colnames(metadata_snmultiome), colnames(metadata_snrnaseq)
)
metadata_snrnaseq$Barcode <- paste0(
  metadata_snrnaseq$Barcode, "_", metadata_snrnaseq$Sample
)
metadata <- rbind(
  metadata_snmultiome[, common_colnames],
  metadata_snrnaseq[, common_colnames]
)

common_cells <- intersect(new_colnames, metadata$Barcode)
log_message(
  "Common cells between object and metadata:", length(common_cells)
)

object <- object[, common_cells]

metadata_filtered <- metadata[metadata$Barcode %in% common_cells, ]
index <- match(colnames(object), metadata_filtered$Barcode)
metadata_filtered <- metadata_filtered[index, ]

object <- AddMetaData(object, metadata_filtered)

age_mapping <- c(
  "HM2013017" = 58,
  "HM20200905" = 44,
  "HM20200927" = 52,
  "HM20201129" = 69,
  "HM20201213" = 47,
  "HM20201222" = 66
)

age_vector <- age_mapping[as.character(object$orig.ident)]
names(age_vector) <- rownames(object@meta.data)
object$age <- age_vector

object$stage <- dplyr::case_when(
  object$age >= 0 & object$age <= 12 ~ "S6", # childhood (0-12 years)
  object$age >= 12 & object$age <= 20 ~ "S7", # adolescence (12-20 years)
  object$age >= 21 & object$age <= 40 ~ "S8", # young adulthood (21-40 years)
  object$age >= 41 & object$age <= 60 ~ "S9", # middle age (41-60 years)
  object$age >= 61 ~ "S10", # old age (61+ years)
  TRUE ~ NA_character_
)
object$dataset <- "PRJCA015229"
object$dataset_number <- "D22"

saveRDS(
  object,
  file.path(
    res_dir, "2024CellGenomics_snRNA-seq.rds"
  )
)
