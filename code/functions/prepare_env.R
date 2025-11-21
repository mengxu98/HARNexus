packages <- c(
  "optparse",
  "Boruta",
  "lightgbm",
  "caret",
  "ggplot2",
  "dplyr",
  "tidyr",
  "stringr",
  "tibble",
  "readr",
  "forcats",
  "purrr",
  "recipes",
  "workflows",
  "tune",
  "yardstick",
  "glmnet",
  "xgboost",
  "tidyverse",
  "biomaRt",
  "rtracklayer",
  "igraph",
  "ggraph",
  "viridis",
  "RColorBrewer",
  "ggrepel",
  "ComplexHeatmap",
  "circlize",
  "UpSetR",
  "data.table",
  "ggplot2",
  "ggsignif",
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot",
  "patchwork",
  "tidyr",
  "igraph",
  "ggraph",
  "uwot",
  "ggwordcloud",
  "caret",
  "glmnet",
  "Metrics",
  "xgboost",
  "lightgbm",
  "ggplot2",
  "reshape2",
  "dplyr",
  "gridExtra",
  "cowplot",
  "patchwork",
  "corrplot",
  "tidyr",
  "rvest",
  "stringr",
  "tidyr",
  "rvest",
  "stringr",
  "openxlsx",
  "Seurat",
  "ggplot2",
  "harmony",
  "mengxu98/lisi",
  "patchwork",
  "ggpubr",
  "lubridate",
  "VennDiagram",
  "ggplot2",
  "ggpointdensity",
  "davidsjoberg/ggsankey",
  "clusterProfiler",
  "org.Hs.eg.db",
  "dplyr",
  "gridExtra",
  "grid",
  "RColorBrewer",
  "patchwork",
  "ggVennDiagram",
  "ggforce",
  "flextable",
  "officer",
  "AUCell",
  "Biostrings",
  "GenomicRanges",
  "TFBSTools",
  "JASPAR2024",
  "motifmatchr",
  "universalmotif",
  "eulerr",
  "BSgenome.Hsapiens.UCSC.hg38",
  "BSgenome.Ptroglodytes.UCSC.panTro5"
)

if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}
if (!requireNamespace("scop", quietly = TRUE)) {
  pak::pak("mengxu98/scop")
}
if (!requireNamespace("thisutils", quietly = TRUE)) {
  pak::pak("thisutils")
}

log_message <- thisutils::log_message

library(scop)

load_packages <- function(
    packages,
    force_update = FALSE) {
  log_message("Installing and loading packages...")
  for (pkg in packages) {
    pkg_name <- if (grepl("/", pkg)) sub(".*/", "", pkg) else pkg

    if (force_update || !requireNamespace(pkg_name, quietly = TRUE)) {
      log_message(paste("Installing", pkg, "..."))
      pak::pak(pkg, upgrade = force_update)
    }
    suppressMessages(
      suppressWarnings(
        suppressPackageStartupMessages(
          library(pkg_name, character.only = TRUE)
        )
      )
    )
  }
  log_message("All packages installed and loaded")
}

load_packages(packages)

check_dir <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    log_message(
      "{.path {dir_path}} does not exist. Creating it"
    )
    dir.create(dir_path, recursive = TRUE)
  }
  return(dir_path)
}
