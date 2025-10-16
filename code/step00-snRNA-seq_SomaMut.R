rm(list = ls())
gc()

# papaer: https://doi.org/10.1038/s41586-025-09435-8
# code: ~
# data: https://publications.wenglab.org/SomaMut/

source("code/functions/packages.R")

data_dir <- "../../data/BrainData/SomaMut"
res_dir <- "../../data/BrainData/results/SomaMut/"
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

log_message("Start loading data...")

object <- readRDS(
  file.path(
    data_dir, "pfc.clean.rds"
  )
)
table(object$age)
object$stage <- dplyr::case_when(
  object$age >= 0 & object$age <= 12 ~ "S6",     # childhood (0-12 years)
  object$age >= 12 & object$age <= 20 ~ "S7",    # adolescence (12-20 years)
  object$age >= 21 & object$age <= 40 ~ "S8",    # young adulthood (21-40 years)
  object$age >= 41 & object$age <= 60 ~ "S9",    # middle age (41-60 years)
  object$age >= 61 ~ "S10",                      # old age (61+ years)
  TRUE ~ NA_character_
)

object$dataset <- "SomaMut"
object$dataset_number <- "D24"

log_message("Save data...")
saveRDS(
  object,
  file.path(res_dir, "SomaMut_snRNA-seq.rds")
)
