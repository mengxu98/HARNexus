packages <- c(
  "AUCell",
  "Biostrings",
  "BSgenome.Hsapiens.UCSC.hg38",
  "BSgenome.Ptroglodytes.UCSC.panTro5",
  "JASPAR2024",
  "Metrics",
  "RColorBrewer",
  "Seurat",
  "TFBSTools",
  "UpSetR",
  "VennDiagram",
  "biomaRt",
  "caret",
  "circlize",
  "clusterProfiler",
  "ComplexHeatmap",
  "cowplot",
  "corrplot",
  "data.table",
  "davidsjoberg/ggsankey",
  "dplyr",
  "enrichplot",
  "eulerr",
  "flextable",
  "forcats",
  "GenomicRanges",
  "ggforce",
  "ggraph",
  "ggpointdensity",
  "ggpubr",
  "ggrepel",
  "ggsignif",
  "ggVennDiagram",
  "ggwordcloud",
  "ggplot2",
  "glmnet",
  "grid",
  "gridExtra",
  "harmony",
  "igraph",
  "inferCSN",
  "lightgbm",
  "lubridate",
  "Matrix",
  "matrixStats",
  "mengxu98/lisi",
  "mengxu98/thisplot",
  "motifmatchr",
  "officer",
  "openxlsx",
  "optparse",
  "org.Hs.eg.db",
  "patchwork",
  "purrr",
  "readr",
  "recipes",
  "reshape2",
  "rtracklayer",
  "rvest",
  "stringr",
  "tidyr",
  "tidyverse",
  "tibble",
  "tune",
  "universalmotif",
  "uwot",
  "viridis",
  "xgboost",
  "yardstick"
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

color_sets <- attr(thisplot::chinese_colors, "color_sets", exact = TRUE)

colors32 <- color_sets$ChineseSet32
colors128 <- color_sets$ChineseSet128
colors2 <- c("gray80", "#15559A")

colors9 <- c(
  "#8076A3",
  "#ED5736",
  "#0AA344",
  "#2177B8",
  "#D70440",
  "#F9BD10",
  "#B14B28",
  "#006D87",
  "#5E7987"
)

color_celltype <- c(
  "Radial glia" = "#8076A3",
  "Neuroblasts" = "#ED5736",
  "Excitatory neurons" = "#0AA344",
  "Inhibitory neurons" = "#2177B8",
  "Astrocytes" = "#D70440",
  "Oligodendrocyte progenitor cells" = "#F9BD10",
  "Oligodendrocytes" = "#B14B28",
  "Microglia" = "#006D87",
  "Endothelial cells" = "#5E7987"
)

marker_genes <- c(
  # Radial glia
  "PAX6", "VIM", "GLI3",
  # Endothelial cells
  "CLDN5", "PECAM1", "VWF", "FLT1",
  # Inhibitory neurons
  "GAD1", "GAD2", "SLC6A1",
  # Oligodendrocyte progenitor cells (OPCs)
  "PDGFRA", "CSPG4", "OLIG1", "OLIG2", "SOX10",
  # Microglia
  "CX3CR1", "P2RY12", "CSF1R",
  # Neuroblasts
  "STMN2",
  # Excitatory neurons
  "SLC17A7", "CAMK2A", "SATB2",
  # Astrocytes
  "GFAP", "AQP4", "ALDH1L1", "FGFR3", "GJA1",
  # Oligodendrocytes
  "MOG", "MAG", "CLDN11"
)
