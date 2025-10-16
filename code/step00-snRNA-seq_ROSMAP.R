rm(list = ls())
gc()

source("code/functions/packages.R")

# paper: https://doi.org/10.1016/j.cell.2023.08.039
# code: https://github.com/mathyslab7/ROSMAP_snRNAseq_PFC/
# data: https://compbio.mit.edu/ad_aging_brain/

system("bash code/functions/download_rosmap_rds-data.sh")
system("bash code/functions/download_rosmap_ucsc_snRNA-seq.sh")
system("bash code/functions/download_rosmap_ucsc_snRNA-seq&snATAC-seq.sh")
system("bash code/functions/download_rosmap_ucsc_snATAC-seq.sh")
system("bash code/functions/download_rosmap_ucsc_snATAC-seq_Epigenomic.sh")


data_dir <- "../../data/BrainData/ROSMAP"
res_dir <- "../../data/BrainData/results/ROSMAP/"

# download metadata:
download.file(
  "https://personal.broadinstitute.org/cboix/ad427_data/Data/Metadata/individual_metadata_deidentified.tsv",
  file.path(data_dir, "individual_metadata_deidentified.tsv")
)

dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

log_message("Start loading data...")

matrix <- Read10X(
  file.path(
    data_dir, "RNA"
  ),
  gene.column = 1
)

meta_data <- read.csv(
  file.path(
    data_dir, "RNA/meta.tsv"
  ),
  sep = "\t"
)
meta_data <- meta_data[meta_data$ADdiag3types == "nonAD", ]
intersect_cells <- intersect(colnames(matrix), meta_data$cellId)
matrix <- matrix[, intersect_cells]
meta_data <- meta_data[meta_data$cellId %in% intersect_cells, ]
object1 <- CreateSeuratObject(
  counts = matrix,
  project = "ROSMAP",
  meta.data = meta_data
)

matrix <- Read10X(
  file.path(
    data_dir, "RNA"
  ),
  gene.column = 1
)

meta_data <- read.csv(
  file.path(
    data_dir, "RNA+ATAC/meta.tsv"
  ),
  sep = "\t"
)
meta_data <- meta_data[meta_data$ADdiag3types == "nonAD", ]
intersect_cells <- intersect(colnames(matrix), meta_data$cellId)
matrix <- matrix[, intersect_cells]
meta_data <- meta_data[meta_data$cellId %in% intersect_cells, ]
object2 <- CreateSeuratObject(
  counts = matrix,
  project = "ROSMAP",
  meta.data = meta_data
)

log_message("Save data...")
saveRDS(
  object,
  file.path(res_dir, "ROSMAP_nonAD_snRNA-seq.rds")
)
