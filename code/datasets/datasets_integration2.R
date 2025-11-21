rm(list = ls())
gc()

source("code/functions/prepare_env.R")

data_dir <- "../../data/BrainData/processed/"
res_dir <- check_dir("../../data/BrainData/integration/")

datasets <- c(
  "AllenM1", "EGAD00001006049", "EGAS00001006537",
  "GSE103723", "GSE104276", "GSE144136", "GSE168408", "GSE178175",
  "GSE186538", "GSE199762", "GSE202210", "GSE204683", "GSE207334",
  "GSE212606", "GSE217511", "GSE67835", "GSE81475", "GSE97942",
  "Li_et_al_2018", "Ma_et_al_2022", "Nowakowski_et_al_2017",
  "PRJCA015229", "ROSMAP", "SomaMut"
)

if (!file.exists(file.path(res_dir, "objects_list.rds"))) {
  objects_list <- list()
  for (dataset in datasets) {
    log_message("Loading: {.val {dataset}}")
    file_path <- file.path(data_dir, dataset, paste0(dataset, "_processed.rds"))
    objects_list[[dataset]] <- readRDS(file_path)
  }
  for (dataset in datasets) {
    object <- objects_list[[dataset]]
    object$unit <- ifelse(grepl("PCW", object$Age), "PCW", "Years")
    object$Stage <- dplyr::case_when(
      object$Age >= 4 & object$Age < 8 & object$unit == "PCW" ~ "Embryonic (4–8 PCW)",
      object$Age >= 8 & object$Age < 10 & object$unit == "PCW" ~ "Early fetal (8–10 PCW)",
      object$Age >= 10 & object$Age < 13 & object$unit == "PCW" ~ "Early fetal (10–13 PCW)",
      object$Age >= 13 & object$Age < 16 & object$unit == "PCW" ~ "Early mid-fetal (13–16 PCW)",
      object$Age >= 16 & object$Age < 19 & object$unit == "PCW" ~ "Early mid-fetal (16–19 PCW)",
      object$Age >= 19 & object$Age < 24 & object$unit == "PCW" ~ "Late mid-fetal (19–24 PCW)",
      object$Age >= 24 & object$Age < 38 & object$unit == "PCW" ~ "Late fetal (24–38 PCW)",
      object$Age >= 0 & object$Age < 0.5 & object$unit == "Years" ~ "Neonatal & early infancy (0–0.5Y)",
      object$Age >= 0.5 & object$Age < 1 & object$unit == "Years" ~ "Late infancy (0.5–1Y)",
      object$Age >= 1 & object$Age < 6 & object$unit == "Years" ~ "Early childhood (1–6Y)",
      object$Age >= 6 & object$Age < 12 & object$unit == "Years" ~ "Middle & late childhood (6–12Y)",
      object$Age >= 12 & object$Age < 20 & object$unit == "Years" ~ "Adolescence (12–20Y)",
      object$Age >= 20 & object$Age < 40 & object$unit == "Years" ~ "Young adulthood (20–40Y)",
      object$Age >= 40 & object$Age < 60 & object$unit == "Years" ~ "Middle adulthood (40–60Y)",
      object$Age >= 60 & object$unit == "Years" ~ "Late adulthood (60+Y)",
      TRUE ~ NA_character_
    )
    objects_list[[dataset]] <- object
  }
  saveRDS(objects_list, file.path(res_dir, "objects_list.rds"))
} else {
  objects_list <- readRDS(file.path(res_dir, "objects_list.rds"))
}

common_genes <- Reduce(intersect, lapply(objects_list, rownames))
for (dataset in datasets) {
  obj <- objects_list[[dataset]]
  obj <- obj[common_genes, ]
  objects_list[[dataset]] <- obj
}

saveRDS(
  objects_list, file.path(res_dir, "objects_list_common_genes.rds")
)
objects <- merge(
  objects_list[[1]],
  y = objects_list[-1]
)
rm(objects_list)
gc()
objects <- JoinLayers(objects)


if (requireNamespace("future", quietly = TRUE)) {
  library(future)
}
options(future.globals.maxSize = Inf)
plan("sequential")

for (dataset in datasets) {
  obj <- objects_list[[dataset]]
  log_message("Processing: {.val {dataset}}")
  obj <- JoinLayers(obj)
  obj[["percent.mt"]] <- PercentageFeatureSet(
    obj,
    pattern = "^MT-|^mt-"
  )
  obj <- subset(obj,
    subset =
      nFeature_RNA > 300 &
        nCount_RNA > 500 &
        percent.mt < 15
  )
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 3000)
  obj <- RunPCA(obj, npcs = 30, approx = TRUE, features = VariableFeatures(obj))

  objects_list[[dataset]] <- obj
}


mini_list <- list()

for (dataset in datasets) {
  obj <- objects_list[[dataset]]

  obj <- FindNeighbors(obj, dims = 1:30, k.param = 10)
  obj <- FindClusters(obj, resolution = 40) # produces mini-clusters

  obj$mini_cluster <- obj$seurat_clusters

  # Average expression per mini-cluster
  avg <- AverageExpression(obj, features = VariableFeatures(obj))$RNA

  mini_list[[dataset]] <- avg
}
mini_mat <- do.call(cbind, mini_list)
saveRDS(mini_mat, file.path(res_dir, "mini_mat.rds"))

library(harmony)

mini_obj <- CreateSeuratObject(mini_mat)
mini_obj <- NormalizeData(mini_obj)
mini_obj <- FindVariableFeatures(mini_obj, nfeatures = 1500)
mini_obj <- ScaleData(mini_obj)
mini_obj <- RunPCA(mini_obj, approx = TRUE)

mini_obj <- RunHarmony(mini_obj, group.by.vars = "dataset")

mini_obj <- RunUMAP(mini_obj, reduction = "harmony", dims = 1:30)
for (dataset in datasets) {
  obj <- objects_list[[dataset]]

  # Map mini-cluster UMAP to each cell
  mini_ids <- obj$seurat_clusters
  obj$UMAP_1 <- mini_obj@reductions$umap@cell.embeddings[mini_ids, 1]
  obj$UMAP_2 <- mini_obj@reductions$umap@cell.embeddings[mini_ids, 2]

  objects_list[[dataset]] <- obj
}




for (dataset in datasets) {
  log_message("Processing: {.val {dataset}}")
  obj <- objects_list[[dataset]]
  obj <- FindNeighbors(obj, dims = 1:30)
  obj <- FindClusters(obj, resolution = 2.0)
  DefaultAssay(obj) <- "SCT"
  markers <- FindAllMarkers(obj, only.pos = TRUE)
  bad_cells <- WhichCells(
    obj,
    expression =
      percent.mt > 20 |
        grepl("^HSP", rownames(obj), ignore.case = TRUE)
  )
  obj <- subset(obj, cells = setdiff(colnames(obj), bad_cells))

  objects_list[[dataset]] <- obj
}

library(matrixStats)

percentile_list <- list()

for (dataset in datasets) {
  obj <- objects_list[[dataset]]

  obj <- SCTransform(obj, vars.to.regress = c("percent.mt"))
  obj <- RunPCA(obj)
  obj <- FindNeighbors(obj, dims = 1:30)
  obj <- FindClusters(obj, resolution = 1.0)

  expr <- GetAssayData(obj, layer = "scale.data")
  v <- rowVars(expr)
  rank <- rank(v) / length(v) * 100
  percentile_list[[dataset]] <- rank
}

percentile_mat <- do.call(cbind, percentile_list)
global_rank <- rowMedians(percentile_mat)

inform_genes <- names(sort(global_rank, decreasing = TRUE))[1:1500]

blacklist <- c(
  grep("^MT-", rownames(percentile_mat), value = TRUE),
  grep("^RPS", rownames(percentile_mat), value = TRUE),
  grep("^RPL", rownames(percentile_mat), value = TRUE),
  "MALAT1"
)
inform_genes <- setdiff(inform_genes, blacklist)

mini_list <- list()
for (dataset in datasets) {
  obj <- objects_list[[dataset]]

  obj_use <- obj[inform_genes, ]
  obj_use <- RunPCA(obj_use)

  obj_use <- FindNeighbors(obj_use, dims = 1:30, k.param = 10)
  obj_use <- FindClusters(obj_use, resolution = 50)

  obj_use$mini <- obj_use$seurat_clusters

  avg_exp <- AverageExpression(
    obj_use,
    assays = "SCT",
    features = inform_genes
  )$SCT
  mini_list[[dataset]] <- avg_exp
}

library(harmony)
mini_obj <- CreateSeuratObject(counts = mini_mat)
mini_obj <- ScaleData(mini_obj)
mini_obj <- RunPCA(mini_obj, features = rownames(mini_mat))

mini_obj <- RunHarmony(mini_obj, group.by.vars = "dataset")
mini_obj <- RunUMAP(mini_obj, dims = 1:50, reduction = "harmony")

for (dataset in datasets) {
  obj <- objects_list[[dataset]]
  mini <- mini_list[[dataset]]

  obj$mini_id <- obj$seurat_clusters

  obj$harmony_umap1 <- mini_obj$umap_1[obj$mini_id]
  obj$harmony_umap2 <- mini_obj$umap_2[obj$mini_id]

  objects_list[[dataset]] <- obj
}
