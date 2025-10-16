rm(list = ls())
gc()

source("code/functions/packages.R")

# paper: https://doi.org/10.1038/s41586-024-08504-8
# code:
#   https://github.com/lsteuernagel/HYPOMAP
#   https://github.com/georgiedowsett/HYPOMAP
#   https://github.com/lsteuernagel/scIntegration
#   https://github.com/mrcepid-rap
# data: https://cellxgene.cziscience.com/collections/d0941303-7ce3-4422-9249-cf31eb98c480
# data(spatial): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278848

data_dir <- "../../data/BrainData/HYPOMAP"
res_dir <- "../../data/BrainData/results/HYPOMAP/"
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

log_message("Start loading data...")

object <- readRDS(
  file.path(
    data_dir, "human_HYPOMAP_snRNASeq.rds"
  )
)
table(object$age_years)

# Add stage column based on age
object$stage <- dplyr::case_when(
  object$age_years >= 0 & object$age_years <= 12 ~ "S6",     # childhood (0-12 years)
  object$age_years >= 12 & object$age_years <= 20 ~ "S7",    # adolescence (12-20 years)
  object$age_years >= 21 & object$age_years <= 40 ~ "S8",    # young adulthood (21-40 years)
  object$age_years >= 41 & object$age_years <= 60 ~ "S9",    # middle age (41-60 years)
  object$age_years >= 61 ~ "S10",                            # old age (61+ years)
  TRUE ~ NA_character_
)
object$dataset <- "HYPOMAP"
object$dataset_number <- "D23"

log_message("Save data...")
saveRDS(
  object,
  file.path(res_dir, "HYPOMAP_snRNA-seq.rds")
)
