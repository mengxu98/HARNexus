rm(list = ls())
gc()

source("code/functions/prepare_env.R")

data_dir <- "../../data/BrainData/raw/"
res_dir <- check_dir("../../data/BrainData/processed/")

log_message("Start loading data...")

datasets <- c(
  "GSE103723", # TPM
  "GSE104276",
  # "GSE126836", # no age information
  "GSE186538", "GSE199762", "GSE204683",
  "GSE212606", "GSE217511",
  "GSE67835", # 416 cells
  # "GSE71585", # mouse data?
  # "GSE76381", # less cells
  "GSE81475",
  "GSE97942", "Li_et_al_2018", "Nowakowski_et_al_2017",
  "SRP041736" # long non-coding RNAs
)

object_list <- list()
not_found <- c()
for (dataset in datasets) {
  log_message("Processing: {.val {dataset}}")
  data_file <- paste0(res_dir, dataset, "/", dataset, ".rds")
  if (!file.exists(data_file)) {
    log_message("File not found: {.val {data_file}}")
    not_found <- c(not_found, dataset)
    next
  }
  seurat <- readRDS(data_file)
  object_list[[dataset]] <- seurat
  rm(seurat)
  gc()
}


if (!file.exists(paste0(res_dir_harnexus, "/merged_seurat_object.rds"))) {
  seurat_objects <- list()

  data_files <- list.files(
    data_dir_harnexus,
    pattern = ".*_seurat\\.Rdata$",
    full.names = TRUE
  )

  all_metadata <- list()
  all_counts <- list()

  for (file in data_files) {
    sample_name <- gsub(".*/(.*?)_seurat\\.Rdata$", "\\1", file)
    load(file)

    if (exists("seurat") && inherits(seurat, "Seurat")) {
      log_message("Processing: {.val {sample_name}}")

      all_counts[[sample_name]] <- GetAssayData(
        seurat,
        layer = "counts"
      )

      meta_data <- seurat@meta.data
      meta_data$sample_id <- sample_name
      all_metadata[[sample_name]] <- meta_data

      rm(seurat)
      gc()
    }
  }

  combined_metadata <- do.call(rbind, all_metadata)

  log_message("Creating merged Seurat object...")
  object_harnexus <- CreateSeuratObject(counts = all_counts)
  object_harnexus@meta.data <- combined_metadata
  rownames(object_harnexus@meta.data) <- Cells(object_harnexus)
  rm(
    all_counts,
    all_metadata,
    combined_metadata
  )
  gc()

  # add dataset information
  dataset_info <- c(
    D1 = "Nowakowski_et_al_2017", D2 = "GSE67835",
    D3 = "GSE104276", D4 = "GSE103723",
    D5 = "Li_et_al_2018", D6 = "Li_et_al_2018",
    D7 = "GSE97942", D8 = "GSE76381",
    D9 = "GSE71585", D10 = "GSE126836",
    D11 = "SRP041736", D12 = "GSE81475",
    D13 = "AllenBrainAtlas", D14 = "GSE212606",
    D15 = "GSE204683", D16 = "GSE199762",
    D17 = "GSE217511", D18 = "GSE186538",
    D19 = "GSE199243", D20 = "GSE168408",
    D21 = "GSE140231"
  )
  brain_dataset_vec <- as.character(object_harnexus$dataset)
  mapped_datasets <- dataset_info[brain_dataset_vec]
  names(mapped_datasets) <- colnames(object_harnexus)
  object_harnexus$Dataset_ID <- mapped_datasets
  na_indices <- is.na(object_harnexus$Dataset_ID)
  object_harnexus$Dataset_ID[na_indices] <- brain_dataset_vec[na_indices]
  filtered_cells <- colnames(
    object_harnexus
  )[object_harnexus$Dataset_ID != "GSE168408"]
  object_harnexus <- object_harnexus[, filtered_cells]
  filtered_cells <- colnames(
    object_harnexus
  )[object_harnexus$Dataset_ID != "AllenBrainAtlas"]
  object_harnexus <- object_harnexus[, filtered_cells]
  filtered_cells <- colnames(
    object_harnexus
  )[object_harnexus$Dataset_ID != "GSE199243"]
  object_harnexus <- object_harnexus[, filtered_cells]
  object_harnexus$Dataset <- object_harnexus$Dataset_ID
  object_harnexus$Region <- object_harnexus$Area
  object_harnexus$CellType <- object_harnexus$cell_name

  log_message("Unifying Age columns by dataset...")

  age_unified <- rep(NA_character_, ncol(object_harnexus))
  datasets <- unique(object_harnexus$Dataset_ID)

  for (dataset_id in datasets) {
    dataset_indices <- object_harnexus$Dataset_ID == dataset_id
    dataset_cells <- colnames(object_harnexus)[dataset_indices]
    dataset_subset <- object_harnexus[, dataset_cells]

    dataset_age <- rep(NA_character_, ncol(dataset_subset))

    if (dataset_id == "GSE67835") {
      age <- dataset_subset@meta.data$Age
      age <- as.character(age)
      week_indices <- !is.na(age) & grepl("16-18w", age)
      dataset_age[week_indices] <- "16-18 PCW"
      other_indices <- !is.na(age) & !week_indices
      dataset_age[other_indices] <- age[other_indices]
    } else if (dataset_id == "Li_et_al_2018") {
      if ("Age_in_Weeks" %in% colnames(dataset_subset@meta.data)) {
        age_weeks <- as.numeric(as.character(dataset_subset@meta.data$Age_in_Weeks))
        week_indices <- !is.na(age_weeks) & age_weeks > 0
        dataset_age[week_indices] <- paste0(age_weeks[week_indices], " PCW")
      }
      if ("Age_in_days" %in% colnames(dataset_subset@meta.data)) {
        age_days <- as.numeric(as.character(dataset_subset@meta.data$Age_in_days))
        day_indices <- !is.na(age_days) & age_days > 0 & is.na(dataset_age)
        dataset_age[day_indices] <- paste0(round(age_days[day_indices] / 365.25, 2), " years")
      }
    } else if (dataset_id %in% c(
      "GSE81475", "GSE103723", "Nowakowski_et_al_2017",
      "SRP041736", "GSE104276", "GSE76381"
    )) {
      if ("Age_in_Weeks" %in% colnames(dataset_subset@meta.data)) {
        age_weeks <- as.numeric(as.character(dataset_subset@meta.data$Age_in_Weeks))
        week_indices <- !is.na(age_weeks) & age_weeks > 0
        dataset_age[week_indices] <- paste0(age_weeks[week_indices], " PCW")
      }
    } else if (dataset_id %in% c("GSE97942", "GSE71585")) {
      if ("Age_in_days" %in% colnames(dataset_subset@meta.data)) {
        age_days <- as.numeric(as.character(dataset_subset@meta.data$Age_in_days))
        day_indices <- !is.na(age_days) & age_days > 0
        dataset_age[day_indices] <- paste0(round(age_days[day_indices] / 365.25, 2), " years")
      }
    } else {
      if ("Age_in_Weeks" %in% colnames(dataset_subset@meta.data)) {
        age_weeks <- as.numeric(as.character(dataset_subset@meta.data$Age_in_Weeks))
        week_indices <- !is.na(age_weeks) & age_weeks > 0
        dataset_age[week_indices] <- paste0(age_weeks[week_indices], " PCW")
      }
      if ("Age" %in% colnames(dataset_subset@meta.data)) {
        age <- as.character(dataset_subset@meta.data$Age)
        age[age == "NA" | age == ""] <- NA
        age_indices <- !is.na(age) & is.na(dataset_age)
        dataset_age[age_indices] <- age[age_indices]
      }
      if ("Age_in_days" %in% colnames(dataset_subset@meta.data)) {
        age_days <- as.numeric(as.character(dataset_subset@meta.data$Age_in_days))
        day_indices <- !is.na(age_days) & age_days > 0 & is.na(dataset_age)
        dataset_age[day_indices] <- paste0(round(age_days[day_indices] / 365.25, 2), " years")
      }
    }

    age_unified[dataset_indices] <- dataset_age
    log_message(
      "Processed {.val {dataset_id}}: {.val {sum(!is.na(dataset_age))}} / {.val {length(dataset_age)}} cells with Age info"
    )

    rm(dataset_subset, dataset_age)
    gc()
  }

  object_harnexus$Age <- age_unified

  log_message("Age columns unified")

  GSE186538 <- subset(object_harnexus, subset = Dataset_ID == "GSE186538")
  object_harnexus <- object_harnexus[, !colnames(object_harnexus) %in% colnames(GSE186538)]
  GSE186538 <- GSE186538[, GSE186538$orig.ident != "D18_HSB628"]
  GSE186538_data <- c(D18_HSB231 = 79, D18_HSB237 = 51)
  GSE186538_meta <- data.frame(
    Age = GSE186538_data[GSE186538$orig.ident],
    Region = "entorhinal cortex",
    stringsAsFactors = FALSE
  )
  rownames(GSE186538_meta) <- colnames(GSE186538)
  GSE186538 <- AddMetaData(GSE186538, GSE186538_meta)
  object_harnexus <- merge(object_harnexus, GSE186538)

  GSE212606 <- subset(object_harnexus, subset = Dataset_ID == "GSE212606")
  object_harnexus <- object_harnexus[, !colnames(object_harnexus) %in% colnames(GSE212606)]
  GSE212606_meta <- read.csv(
    "../../data/BrainData/raw/GSE212606/GSM6657986_cell_annotation.csv",
    header = TRUE,
    row.names = 1
  )
  common_cells <- intersect(colnames(GSE212606), rownames(GSE212606_meta))
  GSE212606_meta <- GSE212606_meta[common_cells, ]
  GSE212606_meta$Cells <- rownames(GSE212606_meta)
  GSE212606_data <- data.frame(
    Individual_ID = c("5459", "5356", "1311", "1306", "1304", "1247"),
    Age = c(70, 94, 85, 83, 81, 94),
    Sex = c("Female", "Male", "Female", "Female", "Male", "Male"),
    stringsAsFactors = FALSE
  )
  GSE212606_meta <- merge(GSE212606_meta, GSE212606_data, by = "Individual_ID")
  rownames(GSE212606_meta) <- GSE212606_meta$Cells
  GSE212606 <- GSE212606[, GSE212606_meta$Cells]
  GSE212606 <- AddMetaData(GSE212606, GSE212606_meta)
  object_harnexus <- merge(object_harnexus, GSE212606)

  GSE140231 <- subset(object_harnexus, subset = Dataset_ID == "GSE140231")
  object_harnexus <- object_harnexus[, !colnames(object_harnexus) %in% colnames(GSE140231)]
  GSE204683 <- subset(object_harnexus, subset = Dataset_ID == "GSE204683")
  object_harnexus <- object_harnexus[, !colnames(object_harnexus) %in% colnames(GSE204683)]
  rm(GSE204683)
  GSE204683_counts <- readRDS(
    "../../data/BrainData/raw/GSE204683/GSE204683_count_matrix.RDS"
  )

  GSE204683_meta <- read.csv(
    "../../data/BrainData/raw/GSE204683/GSE204683_barcodes.tsv",
    header = TRUE,
    row.names = NULL,
    sep = "\t"
  )

  donor_id_mapping <- c(
    "4" = "LaFet1", "8" = "LaFet2", "11" = "EaFet1", "16" = "EaFet2",
    "150656" = "Adult2", "150666" = "Adult1", "4413" = "Inf1", "4422" = "Inf2",
    "5936" = "Adol2", "5977" = "Child2", "6007" = "Adol1", "6032" = "Child1"
  )

  donor_id_reverse_mapping <- names(donor_id_mapping)
  names(donor_id_reverse_mapping) <- donor_id_mapping

  GSE204683_meta$Donor_ID_numeric <- donor_id_reverse_mapping[GSE204683_meta$Donor.ID]
  GSE204683_meta$Cells <- paste0(GSE204683_meta$Donor_ID_numeric, "_", GSE204683_meta$Barcode)

  counts_cells <- colnames(GSE204683_counts)
  GSE204683_meta <- GSE204683_meta[GSE204683_meta$Cells %in% counts_cells, ]

  missing_cells <- setdiff(counts_cells, GSE204683_meta$Cells)
  if (length(missing_cells) > 0) {
    log_message("Found {.val {length(missing_cells)}} cells in counts but not in metadata, creating metadata for them...")
    donor_ids_missing <- sub("_.*", "", missing_cells)
    barcodes_missing <- sub("^[^_]+_", "", missing_cells)
    missing_meta <- data.frame(
      Barcode = barcodes_missing,
      Donor.ID = donor_id_mapping[donor_ids_missing],
      Cell.type = NA_character_,
      UMI.counts = NA_real_,
      gene.count = NA_real_,
      raw.reads.count = NA_real_,
      mapped.reads.count = NA_real_,
      Donor_ID_numeric = donor_ids_missing,
      Cells = missing_cells,
      stringsAsFactors = FALSE,
      row.names = missing_cells
    )
    GSE204683_meta <- rbind(GSE204683_meta, missing_meta)
  }

  rownames(GSE204683_meta) <- GSE204683_meta$Cells
  GSE204683_meta <- GSE204683_meta[counts_cells, ]

  GSE204683_sample_info <- data.frame(
    sample_name = c(
      "Adol1", "Child1", "EaFet2", "Adult2", "Child2", "Adult1",
      "Inf1", "LaFet1", "LaFet2", "Adol2", "Inf2", "EaFet1"
    ),
    brain_region = c(
      "BA 9/46", "BA 9/46", "cortical plate", "BM_9/10/46", "BA 9/46", "BM_9/10/46",
      "BA 9/46", "cortical plate", "cortical plate", "BA 9/46", "BA 9/46", "cortical plate"
    ),
    Region = c(
      "Cerebral cortex", "Cerebral cortex", "Cerebral cortex",
      "Cerebral cortex", "Cerebral cortex", "Cerebral cortex",
      "Cerebral cortex", "Cerebral cortex", "Cerebral cortex",
      "Cerebral cortex", "Cerebral cortex", "Cerebral cortex"
    ),
    Age = c("14", "4", "19 PCW", "39", "6", "20", "1", "23 PCW", "24 PCW", "14", "1", "18 PCW"),
    Sex = c("Male", "Male", "Female", "Male", "Female", "Female", "Female", "Male", "Female", "Female", "Male", "Female"),
    Batch = c("1", "1", "1", "1", "2", "2", "2", "2", "3", "3", "3", "3"),
    stringsAsFactors = FALSE
  )

  GSE204683_meta <- merge(GSE204683_meta, GSE204683_sample_info, by.x = "Donor.ID", by.y = "sample_name", all.x = TRUE)
  rownames(GSE204683_meta) <- GSE204683_meta$Cells
  common_cells <- intersect(rownames(GSE204683_meta), colnames(GSE204683_counts))
  GSE204683_counts <- GSE204683_counts[, GSE204683_meta$Cells]
  GSE204683_object <- CreateSeuratObject(counts = GSE204683_counts)
  GSE204683_object <- AddMetaData(GSE204683_object, GSE204683_meta)
  GSE204683_object$Dataset_ID <- "GSE204683"
  meta_data <- object_harnexus@meta.data
  object_harnexus <- CreateSeuratObject(counts = GetAssayData(object_harnexus, layer = "counts"), meta.data = meta_data)
  object_harnexus <- merge(
    object_harnexus, GSE204683_object[rownames(object_harnexus), ]
  )

  log_message("Checking Age information for each dataset...")
  datasets <- unique(object_harnexus$Dataset_ID)
  age_check <- data.frame(
    Dataset = datasets,
    Total_Cells = 0,
    Has_Age = 0,
    stringsAsFactors = FALSE
  )
  for (i in seq_along(datasets)) {
    dataset_name <- datasets[i]
    dataset_cells <- object_harnexus$Dataset_ID == dataset_name
    age_check$Total_Cells[i] <- sum(dataset_cells)
    age_check$Has_Age[i] <- sum(!is.na(object_harnexus$Age[dataset_cells]))
  }
  print(age_check)

  object_harnexus <- JoinLayers(object_harnexus)
  saveRDS(
    object_harnexus,
    paste0(res_dir_harnexus, "/merged_seurat_object.rds")
  )
} else {
  object_harnexus <- readRDS(
    paste0(res_dir_harnexus, "/merged_seurat_object.rds")
  )
}

# BICCN

# BTSatlas
res_dir_btsatlas <- "../../data/BrainData/processed/BTSatlas/"
object_btsatlas <- readRDS(
  paste0(res_dir_btsatlas, "BTSatlas_snRNA-seq.rds")
)

# GSE207334 + sestanlab2022
res_dir_gse207334 <- "../../data/BrainData/processed/GSE207334/"
object_gse207334 <- readRDS(
  paste0(res_dir_gse207334, "GSE207334_snRNA-seq.rds")
)
object_sestanlab2022 <- readRDS(
  paste0(res_dir_gse207334, "sestanlab2022_snRNA-seq.rds")
)

# HYPOMAP
res_dir_hypomap <- "../../data/BrainData/processed/HYPOMAP/"
object_hypomap <- readRDS(
  paste0(res_dir_hypomap, "HYPOMAP_snRNA-seq.rds")
)

# SomaMut 367317 cells
res_dir_soma_mut <- "../../data/BrainData/processed/SomaMut/"
object_soma_mut <- readRDS(
  paste0(res_dir_soma_mut, "SomaMut_snRNA-seq.rds")
)

# rosamp 172683 cells
res_dir_rosmap <- "../../data/BrainData/processed/ROSMAP/"
object_rosmap <- readRDS(
  paste0(res_dir_rosmap, "ROSMAP_nonAD_snRNA-seq.rds")
)

if (!file.exists(paste0(res_dir_harnexus, "/merged_seurat_object_processed.rds"))) {
  object_harnexus <- JoinLayers(object_harnexus)
  object_harnexus <- NormalizeData(object_harnexus)
  object_harnexus <- FindVariableFeatures(object_harnexus)
  object_harnexus <- ScaleData(object_harnexus)
  object_harnexus <- RunPCA(object_harnexus)
  object_harnexus <- FindNeighbors(object_harnexus, dims = 1:50)
  object_harnexus <- FindClusters(object_harnexus, resolution = 2)
  object_harnexus <- RunUMAP(object_harnexus, dims = 1:30)

  object_harnexus <- RunHarmony(
    object_harnexus,
    group.by.vars = "dataset",
    reduction = "pca",
    dims.use = 1:30,
    reduction.save = "harmony"
  )

  object_harnexus <- RunUMAP(
    object_harnexus,
    reduction = "harmony",
    dims = 1:30,
    reduction.name = "umap.harmony"
  )

  object_harnexus <- FindNeighbors(
    object_harnexus,
    reduction = "integrated.harmony",
    dims = 1:30
  )
  object_harnexus <- FindNeighbors(
    object_harnexus,
    reduction = "harmony",
    dims = 1:30
  )

  saveRDS(
    object_harnexus,
    paste0(
      res_dir_harnexus, "/merged_seurat_object_processed.rds"
    )
  )

  pca_raw <- object_harnexus@reductions$pca@cell.embeddings
  pca_harmony <- object_harnexus@reductions$harmony@cell.embeddings
  umap_raw <- object_harnexus@reductions$umap@cell.embeddings
  umap_harmony <- object_harnexus@reductions$umap.harmony@cell.embeddings
  meta_data <- object_harnexus@meta.data
  lisi_data <- list(
    pca_raw,
    pca_harmony,
    umap_raw,
    umap_harmony,
    meta_data
  )
  saveRDS(
    lisi_data,
    paste0(res_dir_harnexus, "/lisi_data.rds")
  )
} else {
  object_harnexus <- readRDS(
    paste0(res_dir_harnexus, "/merged_seurat_object_processed.rds")
  )
}
