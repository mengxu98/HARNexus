files <- list.files(
  "data/grn",
  pattern = "*.Rdata",
  full.names = TRUE
)

science_har <- read.csv(
  "data/genome/science_HAR.csv"
)
head(science_har)
tfs <- read.table("data/regulators.txt")[, 1]
tfdb_result <- read.delim(
  "data/AnimalTFDB_result.txt",
  stringsAsFactors = FALSE
)
tfdb_result <- tfdb_result[tfdb_result$TF %in% tfs, ]
head(tfdb_result)
unique(tfdb_result$TF)

get_tf_har_combinations <- function() {
  tf_har_combinations <- data.frame()
  tfs_used_list <- c()

  for (i in 1:nrow(science_har)) {
    tf <- science_har[i, ]
    tf_chrom <- paste0(
      tf$chrom, ":",
      tf$start, "-",
      tf$end
    )
    har_name <- tf$har
    har_name <- gsub("^ZOO", "", har_name)

    matching_har <- tfdb_result[tfdb_result$Sequence == tf_chrom, ]
    tf_har_combinations <- rbind(
      tf_har_combinations,
      data.frame(
        Chrom = tf$chrom,
        Start = tf$start,
        End = tf$end,
        HAR = har_name,
        TF = unique(matching_har$TF)
      )
    )
  }

  return(tf_har_combinations)
}

tf_har_df <- get_tf_har_combinations()
head(tf_har_df)
nrow(tf_har_df)

length(unique(tf_har_df$TF))
length(unique(tf_har_df$HAR))

write.csv(
  tf_har_df,
  "data/networks/csv/tf_har_combinations_with_chrom.csv",
  row.names = FALSE
)

write.csv(
  tf_har_df[, -c(1:3)],
  "data/networks/csv/tf_har_combinations.csv",
  row.names = FALSE
)


get_components <- function(filename) {
  basename <- gsub(".Rdata$", "", basename(filename))
  parts <- strsplit(basename, "_")[[1]]
  brain_region <- parts[1]
  stage <- gsub("grn.*$", "", parts[2])
  return(list(brain_region = brain_region, stage = stage))
}

result <- list()
for (file in files) {
  load(file)

  components <- get_components(file)
  brain_region <- components$brain_region
  stage <- components$stage

  if (is.null(result[[stage]])) {
    result[[stage]] <- list()
  }

  if (is.null(result[[stage]][[brain_region]])) {
    result[[stage]][[brain_region]] <- list()
  }

  result[[stage]][[brain_region]] <- grn_list
}

saveRDS(result, file = "data/networks/reorganized_networks.rds")
result <- readRDS("data/networks/reorganized_networks.rds")

dir.create("data/networks/csv", recursive = TRUE, showWarnings = FALSE)

convert_network <- function(network_data) {
  all_data <- data.frame()

  for (stage in names(network_data)) {
    stage_data <- network_data[[stage]]

    for (region in names(stage_data)) {
      region_data <- stage_data[[region]]

      for (cell_type in names(region_data)) {
        log_message(
          "Processing ", stage, "-", region, "-", cell_type
        )
        network <- region_data[[cell_type]]
        if (!is.null(network) && nrow(network) > 0) {
          network_df <- data.frame(
            Stage = stage,
            Region = region,
            CellType = cell_type,
            TF = network$regulator,
            Target = network$target,
            Weight = network$weight
          )

          all_data <- rbind(all_data, network_df)

          print(head(all_data))
        }
      }
    }
  }
  return(all_data)
}

network_df <- convert_network(result)

write.csv(
  network_df,
  "data/networks/csv/network_data.csv",
  row.names = FALSE
)
