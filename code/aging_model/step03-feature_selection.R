source("code/functions/utils.R")
source("code/functions/packages.R")
source("code/functions/aging_data.R")
source("code/functions/aging_model.R")
source("code/functions/aging_plots.R")

option_list <- list(
  make_option(c("-f", "--feature_file"),
    type = "character",
    default = "results/networks_pfc/s8_s9_common_targets.csv",
    help = "Path to the input feature list CSV file [default %default]",
    metavar = "character"
  ),
  make_option(c("-d", "--data_file"),
    type = "character",
    default = "data/ukb/olink_data_imputed_mean_drop30.csv",
    help = "Path to the main data CSV file [default %default]",
    metavar = "character"
  ),
  make_option(c("-r", "--response_file"),
    type = "character",
    default = "data/ukb/granular_age_april_07_2025.csv",
    help = "Path to the response variable CSV file [default %default]",
    metavar = "character"
  ),
  make_option(c("-o", "--output_dir"),
    type = "character",
    default = "results/aging_model/features_selection/har/",
    help = "Base directory for results [default %default]",
    metavar = "character"
  ),
  make_option(c("-t", "--train_ratio"),
    type = "double",
    default = 0.7,
    help = "Proportion of data for training set [default %default]",
    metavar = "number"
  ),
  make_option(c("-s", "--seed"),
    type = "integer",
    default = 2025,
    help = "Random seed for reproducibility [default %default]",
    metavar = "integer"
  ),
  make_option(c("-i", "--nrounds"),
    type = "integer",
    default = 1000,
    help = "Maximum LightGBM iterations [default %default]",
    metavar = "integer"
  ),
  make_option(c("-b", "--n_boruta"),
    type = "integer",
    default = 200,
    help = "Number of Boruta iterations [default %default]",
    metavar = "integer"
  ),
  make_option(c("-p", "--features_num_boruta"),
    type = "integer",
    default = 200,
    help = "Number of features to pre-select before Boruta [default %default]",
    metavar = "integer"
  ),
  make_option(c("-u", "--features_num_optimized"),
    type = "integer",
    default = 20,
    help = "Number of features to optimize [default %default]",
    metavar = "integer"
  ),
  make_option(c("-k", "--n_folds"),
    type = "integer",
    default = 5,
    help = "Number of folds for cross-validation [default %default]",
    metavar = "integer"
  )
)

log_message("Starting feature selection process...")

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

data_file <- opt$data_file
response_file <- opt$response_file
feature_file <- opt$feature_file

set.seed(opt$seed)
nrounds <- opt$nrounds
n_boruta <- opt$n_boruta
n_folds <- opt$n_folds
train_ratio <- opt$train_ratio
features_num_boruta <- opt$features_num_boruta
features_num_optimized <- opt$features_num_optimized

output_dir <- check_dir(opt$output_dir)
model_dir <- check_dir(paste0(output_dir, "/models"))
fig_dir <- check_dir(paste0(output_dir, "/plots"))

max_depth_pretrain <- 6
max_depth_boruta <- 6
max_depth_optimized <- 5

params_pretrain <- list(
  objective = "regression",
  metric = c("rmse", "mae", "mse", "r2"),
  learning_rate = 0.05,
  max_depth = max_depth_pretrain,
  num_leaves = 2^(max_depth_pretrain - 1),
  feature_fraction = 0.7,
  bagging_fraction = 0.7,
  verbosity = -1,
  num_threads = 6
)

params_boruta <- list(
  objective = "regression",
  metric = c("rmse", "mae", "mse", "r2"),
  learning_rate = 0.01,
  max_depth = max_depth_boruta,
  num_leaves = 2^(max_depth_boruta - 1),
  feature_fraction = 0.7,
  bagging_fraction = 0.7,
  min_data_in_leaf = 5,
  lambda_l1 = 0.01,
  lambda_l2 = 0.05,
  min_gain_to_split = 0.01,
  verbosity = -1,
  num_threads = 6
)

params_optimized <- list(
  objective = "regression",
  metric = c("rmse", "mae", "mse", "r2"),
  learning_rate = 0.01,
  max_depth = max_depth_optimized,
  num_leaves = 2^(max_depth_optimized - 1),
  feature_fraction = 0.65,
  bagging_fraction = 0.65,
  lambda_l1 = 0.05,
  lambda_l2 = 0.1,
  verbosity = -1,
  num_threads = 6
)

log_message("1. Loading and preparing data...")
data <- prepare_data_new(
  data_file = data_file,
  feature_file = feature_file,
  response_file = response_file,
  preprocessing = FALSE,
  precentage = train_ratio,
  age_at = 0
)

data_sets <- create_train_test_sets(
  data$X,
  data$response_variable,
  data$gender_numeric,
  data$train_index,
  add_covariate = FALSE,
  scale = TRUE
)

log_message(
  sprintf(
    "Loaded %d features, %d samples",
    ncol(data$X), nrow(data$X)
  )
)
log_message(
  sprintf(
    "Training set: %d samples, Test set: %d samples",
    nrow(data_sets$X_train_scaled),
    nrow(data_sets$X_test_scaled)
  )
)

log_message("2. Preliminary feature selection...")
dtrain <- lightgbm::lgb.Dataset(
  data = as.matrix(data_sets$X_train_scaled),
  label = data_sets$y_train
)

pre_model <- lightgbm::lgb.train(
  params = params_pretrain,
  data = dtrain,
  nrounds = nrounds
)

train_preds_pretrain <- predict(
  pre_model,
  as.matrix(data_sets$X_train_scaled)
)
test_preds_pretrain <- predict(
  pre_model,
  as.matrix(data_sets$X_test_scaled)
)

metrics_train_pretrain <- calculate_metrics(
  data_sets$y_train,
  train_preds_pretrain
)
metrics_test_pretrain <- calculate_metrics(
  data_sets$y_test,
  test_preds_pretrain
)

train_r2_pretrain <- metrics_train_pretrain$r2
test_r2_pretrain <- metrics_test_pretrain$r2
train_rmse_pretrain <- metrics_train_pretrain$rmse
test_rmse_pretrain <- metrics_test_pretrain$rmse
train_mse_pretrain <- metrics_train_pretrain$mse
test_mse_pretrain <- metrics_test_pretrain$mse
train_mae_pretrain <- metrics_train_pretrain$mae
test_mae_pretrain <- metrics_test_pretrain$mae
train_correlation_pretrain <- metrics_train_pretrain$correlation
test_correlation_pretrain <- metrics_test_pretrain$correlation

log_message(
  sprintf(
    "Pre-train set: R² = %.4f, RMSE = %.4f, MAE = %.4f, MSE = %.4f, Correlation = %.4f",
    train_r2_pretrain,
    train_rmse_pretrain,
    train_mae_pretrain,
    train_mse_pretrain,
    train_correlation_pretrain
  )
)
log_message(
  sprintf(
    "Pre-test set: R² = %.4f, RMSE = %.4f, MAE = %.4f, MSE = %.4f, Correlation = %.4f",
    test_r2_pretrain,
    test_rmse_pretrain,
    test_mae_pretrain,
    test_mse_pretrain,
    test_correlation_pretrain
  )
)
metrics_pretrain <- data.frame(
  metric = c("R²", "RMSE", "MAE", "MSE", "Correlation"),
  train = c(
    train_r2_pretrain,
    train_rmse_pretrain,
    train_mae_pretrain,
    train_mse_pretrain,
    train_correlation_pretrain
  ),
  test = c(
    test_r2_pretrain,
    test_rmse_pretrain,
    test_mae_pretrain,
    test_mse_pretrain,
    test_correlation_pretrain
  )
)
write.csv(
  metrics_pretrain,
  paste0(output_dir, "/metrics_pretrain.csv"),
  row.names = FALSE
)


importance_pretrain <- lightgbm::lgb.importance(
  pre_model,
  percentage = TRUE
)
write.csv(
  importance_pretrain,
  paste0(output_dir, "/features_importance_pretrain.csv"),
  row.names = FALSE
)
features_pretrain <- importance_pretrain$Feature[seq_len(
  min(features_num_boruta, nrow(importance_pretrain))
)]
write.csv(
  features_pretrain,
  paste0(output_dir, "/features_pretrain.csv"),
  row.names = FALSE
)

x_train_filtered <- data_sets$X_train_scaled[, features_pretrain]
x_test_filtered <- data_sets$X_test_scaled[, features_pretrain]
log_message(
  sprintf(
    "Selected top %d important features for further analysis",
    length(features_pretrain)
  )
)

log_message("3. Executing Boruta feature selection...")

start_time <- Sys.time()

importance_wrapper <- function(matrix, y) {
  dtrain <- lightgbm::lgb.Dataset(
    data = as.matrix(matrix),
    label = y
  )
  model <- lightgbm::lgb.train(
    params = params_boruta,
    data = dtrain,
    nrounds = nrounds
  )

  imp <- lightgbm::lgb.importance(model, percentage = TRUE)
  imp_values <- rep(0, ncol(matrix))
  for (i in seq_len(nrow(imp))) {
    feat_idx <- which(colnames(matrix) == imp$Feature[i])
    if (length(feat_idx) > 0) {
      imp_values[feat_idx] <- imp$Gain[i]
    }
  }
  return(imp_values)
}

boruta_result <- Boruta::Boruta(
  x_train_filtered,
  data_sets$y_train,
  doTrace = 1,
  maxRuns = n_boruta,
  getImp = importance_wrapper
)

confirmed_features <- names(
  boruta_result$finalDecision[boruta_result$finalDecision == "Confirmed"]
)
tentative_features <- names(
  boruta_result$finalDecision[boruta_result$finalDecision == "Tentative"]
)
features_optimized <- c(confirmed_features, tentative_features)
if (length(features_optimized) >= features_num_optimized) {
  features_optimized <- head(features_optimized, features_num_optimized)
}

end_time <- Sys.time()
log_message(
  sprintf(
    "Boruta completed, time: %.2f minutes",
    as.numeric(difftime(end_time, start_time, units = "mins"))
  )
)
log_message(
  sprintf(
    "Selected features: %d (Confirmed: %d, Tentative: %d)",
    length(features_optimized),
    length(confirmed_features),
    length(tentative_features)
  )
)

write.csv(
  features_optimized,
  paste0(output_dir, "/features_optimized.csv"),
  row.names = FALSE
)

log_message("4. Training final model...")
if (length(features_optimized) > 0) {
  x_train_final <- x_train_filtered[, features_optimized, drop = FALSE]
  x_test_final <- x_test_filtered[, features_optimized, drop = FALSE]

  dtrain_final <- lightgbm::lgb.Dataset(
    data = as.matrix(x_train_final),
    label = data_sets$y_train
  )
  dtest_final <- lightgbm::lgb.Dataset(
    data = as.matrix(x_test_final),
    label = data_sets$y_test,
    reference = dtrain_final
  )

  model_optimized <- lightgbm::lgb.train(
    params = params_optimized,
    data = dtrain_final,
    nrounds = nrounds,
    valids = list(test = dtest_final),
    early_stopping_rounds = 50,
    verbose = 0
  )

  saveRDS(model_optimized, paste0(model_dir, "/lightgbm_model.rds"))

  train_preds <- predict(model_optimized, as.matrix(x_train_final))
  test_preds <- predict(model_optimized, as.matrix(x_test_final))

  metrics_train <- calculate_metrics(data_sets$y_train, train_preds)
  metrics_test <- calculate_metrics(data_sets$y_test, test_preds)

  train_r2 <- metrics_train$r2
  test_r2 <- metrics_test$r2
  train_rmse <- metrics_train$rmse
  test_rmse <- metrics_test$rmse
  train_mse <- metrics_train$mse
  test_mse <- metrics_test$mse
  train_mae <- metrics_train$mae
  test_mae <- metrics_test$mae
  train_correlation <- metrics_train$correlation
  test_correlation <- metrics_test$correlation

  log_message(
    sprintf(
      "Training set: R² = %.4f, RMSE = %.4f, MAE = %.4f, MSE = %.4f, Correlation = %.4f",
      train_r2,
      train_rmse,
      train_mae,
      train_mse,
      train_correlation
    )
  )
  log_message(
    sprintf(
      "Test set: R² = %.4f, RMSE = %.4f, MAE = %.4f, MSE = %.4f, Correlation = %.4f",
      test_r2,
      test_rmse,
      test_mae,
      test_mse,
      test_correlation
    )
  )

  log_message("5. Calculating predicted age and age gap...")
  all_x <- rbind(
    data_sets$X_train_scaled[, features_optimized, drop = FALSE],
    data_sets$X_test_scaled[, features_optimized, drop = FALSE]
  )
  all_y <- c(data_sets$y_train, data_sets$y_test)
  age_preds <- predict(model_optimized, as.matrix(all_x))

  age_df <- data.frame(
    eid = rownames(all_x),
    Age = all_y,
    Predict_age = age_preds
  )

  age_df$age_gap <- age_df$Predict_age - age_df$Age

  write.csv(
    age_df,
    paste0(output_dir, "/age_results.csv"),
    row.names = FALSE
  )

  log_message("6. Creating visualizations...")

  age_df_sampled <- age_df %>%
    sample_n(size = 5000)
  age_scatter <- inferCSN::plot_scatter(
    age_df_sampled[, c("Age", "Predict_age")]
  ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      color = "#a02a2a",
      linetype = "dashed"
    ) +
    geom_smooth(
      method = "lm",
      color = "white",
      fill = "white",
      alpha = 0.2
    ) +
    ggtitle(
      paste0(
        "Test set R² = ",
        round(test_r2, 3)
      )
    ) + theme_bw()

  ggsave(
    paste0(fig_dir, "/age_scatter.pdf"),
    age_scatter,
    width = 5,
    height = 4
  )

  importance_optimized <- lightgbm::lgb.importance(
    model = model_optimized,
    percentage = TRUE
  )

  if (nrow(importance_optimized) > 0) {
    plot_data <- head(importance_optimized, 20)
    plot_data$type <- "Boruta"
    plot_data <- plot_data[, c("Feature", "type", "Gain")]
    colnames(plot_data) <- c("regulator", "type", "weight")
    importance_plot <- inferCSN::plot_coefficient(
      plot_data,
      show_values = FALSE
    ) +
      xlab("Feature importance") +
      theme_bw() +
      theme(
        legend.position = "none"
      ) +
      ylim(0, max(plot_data$weight) * 1.2) +
      geom_text(
        aes(label = sprintf("%.4f", weight)),
        hjust = 0,
        vjust = 0.5,
        size = 3,
        color = "black"
      )
    ggsave(
      paste0(fig_dir, "/feature_importance.pdf"),
      importance_plot,
      width = 3,
      height = 4
    )
  }

  metrics_optimized <- data.frame(
    metric = c("R²", "RMSE", "MAE", "MSE", "Correlation"),
    train = c(
      train_r2,
      train_rmse,
      train_mae,
      train_mse,
      train_correlation
    ),
    test = c(
      test_r2,
      test_rmse,
      test_mae,
      test_mse,
      test_correlation
    )
  )
  write.csv(
    metrics_optimized,
    paste0(output_dir, "/metrics_optimized.csv"),
    row.names = FALSE
  )

  write.csv(
    importance_optimized,
    paste0(output_dir, "/features_importance_optimized.csv"),
    row.names = FALSE
  )

  log_message(
    "========================Analysis completed!========================",
    message_type = "success"
  )
} else {
  log_message(
    "Warning: Boruta did not select any features, cannot train final model.",
    message_type = "warning"
  )
}
