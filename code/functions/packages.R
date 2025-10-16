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
  "rtracklayer",
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
  "AUCell"
)
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}
pak::pak("mengxu98/scop")
pak::pak("thisutils")
library(scop)
PrepareEnv()
log_message <- thisutils::log_message

install_and_load_packages <- function(
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
      suppressPackageStartupMessages(
        library(pkg_name, character.only = TRUE)
      )
    )
  }
  log_message("All packages installed and loaded")
}

install_and_load_packages(packages)
