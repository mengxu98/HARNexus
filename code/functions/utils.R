check_dir <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    inferCSN::log_message(
      sprintf(
        "Directory %s does not exist. Creating it.", dir_path
      )
    )
    dir.create(dir_path, recursive = TRUE)
  }
  return(dir_path)
}

calculate_metrics <- function(actual, predicted) {
  pred_vector <- as.vector(predicted)
  mse_val <- Metrics::mse(actual, pred_vector)
  rmse_val <- sqrt(mse_val)
  mae_val <- Metrics::mae(actual, pred_vector)
  r2_val <- 1 - sum((actual - pred_vector)^2) / sum((actual - mean(actual))^2)

  return(
    list(
      correlation = cor(actual, pred_vector),
      mse = mse_val,
      mae = mae_val,
      rmse = rmse_val,
      r2 = r2_val
    )
  )
}
