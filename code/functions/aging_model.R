source("code/functions/utils.R")


train_evaluate_models <- function(
    X_train_scaled,
    X_test_scaled,
    y_train,
    y_test,
    model_name,
    params_list = list(
      lasso = list(
        alpha = 1,
        lambda = c(
          1e-15, 1e-10, 1e-8, 1e-5, 1e-4,
          1e-3, 1e-2, 1, 5, 10, 50
        )
      ),
      elastic_net = list(
        l1_ratios = c(
          0.1, 0.5, 0.7, 0.9, 0.95, 0.99
        ),
        alpha_values = c(
          1e-15, 1e-10, 1e-8, 1e-5, 1e-4,
          1e-3, 1e-2, 1, 5, 10, 50
        )
      ),
      lightgbm = list(
        objective = "regression",
        metric = c("rmse", "mae", "mse", "r2"),
        learning_rate = 0.01,
        max_depth = 5,
        num_leaves = 2^(5 - 1),
        feature_fraction = 0.65,
        bagging_fraction = 0.65,
        lambda_l1 = 0.05,
        lambda_l2 = 0.1,
        verbosity = -1,
        nrounds = 1000,
        early_stopping_rounds = 100,
        num_threads = 6
      ),
      xgboost = list(
        objective = "reg:squarederror",
        eta = 0.01,
        max_depth = 8,
        subsample = 0.8,
        colsample_bytree = 0.8,
        min_child_weight = 3,
        gamma = 1,
        lambda = 1,
        alpha = 0.5,
        nrounds = 1000,
        early_stopping_rounds = 100,
        num_threads = 6
      )
    ),
    nfolds = 5,
    model_use = c("LASSO", "Elastic Net", "LightGBM", "XGBoost"),
    model_dir = "") {
  results <- data.frame(
    model = character(),
    method = character(),
    mse = numeric(),
    rmse = numeric(),
    mae = numeric(),
    r2 = numeric(),
    correlation = numeric(),
    mse_train = numeric(),
    rmse_train = numeric(),
    mae_train = numeric(),
    r2_train = numeric(),
    correlation_train = numeric(),
    stringsAsFactors = FALSE
  )
  dir_path <- check_dir(paste0(model_dir, "/", model_name))

  predictions <- list()
  train_predictions <- list()

  log_message("\nTraining LASSO model...\n")
  lasso_model <- glmnet::cv.glmnet(
    as.matrix(X_train_scaled), y_train,
    alpha = params_list$lasso$alpha,
    lambda = params_list$lasso$lambda,
    nfolds = nfolds
  )

  lambda_lasso <- lasso_model$lambda.min
  lasso_fit <- glmnet::glmnet(
    as.matrix(X_train_scaled), y_train,
    alpha = params_list$lasso$alpha,
    lambda = lambda_lasso,
    trace.it = 1
  )
  saveRDS(
    lasso_fit,
    paste(
      model_dir,
      model_name,
      "/lasso_model.rds",
      sep = ""
    )
  )
  lasso_pred <- predict(
    lasso_fit,
    newx = as.matrix(X_test_scaled)
  )
  lasso_train_pred <- predict(
    lasso_fit,
    newx = as.matrix(X_train_scaled)
  )
  predictions[["LASSO"]] <- lasso_pred
  train_predictions[["LASSO"]] <- lasso_train_pred
  lasso_metrics <- calculate_metrics(y_test, lasso_pred)
  lasso_metrics_train <- calculate_metrics(y_train, lasso_train_pred)

  results <- rbind(
    results,
    data.frame(
      model = model_name,
      method = "LASSO",
      mse = lasso_metrics$mse,
      rmse = lasso_metrics$rmse,
      mae = lasso_metrics$mae,
      r2 = lasso_metrics$r2,
      correlation = lasso_metrics$correlation,
      mse_train = lasso_metrics_train$mse,
      rmse_train = lasso_metrics_train$rmse,
      mae_train = lasso_metrics_train$mae,
      r2_train = lasso_metrics_train$r2,
      correlation_train = lasso_metrics_train$correlation
    )
  )

  log_message("Training Elastic Net model...\n")
  elastic_net_performance <- data.frame(
    alpha = numeric(),
    l1_ratio = numeric(),
    mse = numeric()
  )

  require(doMC)
  registerDoMC(cores = 4)
  for (l1_ratio in params_list$elastic_net$l1_ratios) {
    model <- glmnet::cv.glmnet(
      as.matrix(X_train_scaled),
      y_train,
      alpha = l1_ratio,
      lambda = params_list$elastic_net$alpha_values,
      nfolds = nfolds,
      parallel = TRUE,
      type.measure = "mse"
    )
    lambda_min <- model$lambda.min
    pred <- predict(model, s = lambda_min, newx = as.matrix(X_test_scaled))
    mse_val <- mse(y_test, pred)
    elastic_net_performance <- rbind(
      elastic_net_performance,
      data.frame(alpha = lambda_min, l1_ratio = l1_ratio, mse = mse_val)
    )
  }

  best_elastic_net <- elastic_net_performance[which.min(elastic_net_performance$mse), ]
  best_l1_ratio <- best_elastic_net$l1_ratio
  best_alpha <- best_elastic_net$alpha

  elastic_net_model <- glmnet::cv.glmnet(
    as.matrix(X_train_scaled),
    y_train,
    alpha = best_l1_ratio,
    lambda = params_list$elastic_net$alpha_values,
    nfolds = nfolds,
    parallel = TRUE,
    type.measure = "mse"
  )

  saveRDS(
    elastic_net_model,
    paste(
      model_dir,
      model_name,
      "/elastic_net_model.rds",
      sep = ""
    )
  )
  elastic_net_pred <- predict(
    elastic_net_model,
    s = best_alpha,
    newx = as.matrix(X_test_scaled)
  )
  elastic_net_train_pred <- predict(
    elastic_net_model,
    s = best_alpha,
    newx = as.matrix(X_train_scaled)
  )
  predictions[["Elastic Net"]] <- elastic_net_pred
  train_predictions[["Elastic Net"]] <- elastic_net_train_pred
  elastic_net_metrics <- calculate_metrics(y_test, elastic_net_pred)
  elastic_net_metrics_train <- calculate_metrics(y_train, elastic_net_train_pred)
  results <- rbind(
    results,
    data.frame(
      model = model_name,
      method = "Elastic Net",
      mse = elastic_net_metrics$mse,
      rmse = elastic_net_metrics$rmse,
      mae = elastic_net_metrics$mae,
      r2 = elastic_net_metrics$r2,
      correlation = elastic_net_metrics$correlation,
      mse_train = elastic_net_metrics_train$mse,
      rmse_train = elastic_net_metrics_train$rmse,
      mae_train = elastic_net_metrics_train$mae,
      r2_train = elastic_net_metrics_train$r2,
      correlation_train = elastic_net_metrics_train$correlation
    )
  )

  log_message("Training LightGBM model...\n")
  dtrain <- lightgbm::lgb.Dataset(
    data = as.matrix(X_train_scaled),
    label = y_train
  )
  lgb_model <- lightgbm::lgb.train(
    params = params_list$lightgbm,
    data = dtrain,
    nrounds = params_list$lightgbm$nrounds,
    verbose = 0,
    valids = list(
      valid = lightgbm::lgb.Dataset(
        data = as.matrix(X_test_scaled),
        label = y_test
      )
    ),
    early_stopping_rounds = params_list$lightgbm$early_stopping_rounds
  )
  saveRDS(
    lgb_model,
    paste(
      model_dir,
      model_name,
      "/lightgbm_model.rds",
      sep = ""
    )
  )

  lgb_pred <- predict(lgb_model, newdata = as.matrix(X_test_scaled))
  lgb_train_pred <- predict(
    lgb_model,
    newdata = as.matrix(X_train_scaled)
  )
  predictions[["LightGBM"]] <- lgb_pred
  train_predictions[["LightGBM"]] <- lgb_train_pred
  lgb_metrics <- calculate_metrics(y_test, lgb_pred)
  lgb_metrics_train <- calculate_metrics(y_train, lgb_train_pred)
  results <- rbind(
    results, data.frame(
      model = model_name,
      method = "LightGBM",
      mse = lgb_metrics$mse,
      rmse = lgb_metrics$rmse,
      mae = lgb_metrics$mae,
      r2 = lgb_metrics$r2,
      correlation = lgb_metrics$correlation,
      mse_train = lgb_metrics_train$mse,
      rmse_train = lgb_metrics_train$rmse,
      mae_train = lgb_metrics_train$mae,
      r2_train = lgb_metrics_train$r2,
      correlation_train = lgb_metrics_train$correlation
    )
  )

  log_message("Training XGBoost model...\n")
  dtrain_xgb <- xgb.DMatrix(
    data = as.matrix(X_train_scaled),
    label = y_train
  )
  dtest_xgb <- xgb.DMatrix(
    data = as.matrix(X_test_scaled),
    label = y_test
  )

  xgb_model <- xgb.train(
    params = params_list$xgboost,
    data = dtrain_xgb,
    nrounds = params_list$xgboost$nrounds,
    watchlist = list(eval = dtest_xgb, train = dtrain_xgb),
    early_stopping_rounds = params_list$xgboost$early_stopping_rounds,
    print_every_n = 100,
    verbose = 1
  )

  saveRDS(
    xgb_model,
    paste(
      model_dir,
      model_name,
      "/xgboost_model.rds",
      sep = ""
    )
  )

  xgb_pred <- predict(xgb_model, dtest_xgb)
  xgb_train_pred <- predict(xgb_model, dtrain_xgb)
  predictions[["XGBoost"]] <- xgb_pred
  train_predictions[["XGBoost"]] <- xgb_train_pred
  xgb_metrics <- calculate_metrics(y_test, xgb_pred)
  xgb_metrics_train <- calculate_metrics(y_train, xgb_train_pred)
  results <- rbind(
    results, data.frame(
      model = model_name,
      method = "XGBoost",
      mse = xgb_metrics$mse,
      rmse = xgb_metrics$rmse,
      mae = xgb_metrics$mae,
      r2 = xgb_metrics$r2,
      correlation = xgb_metrics$correlation,
      mse_train = xgb_metrics_train$mse,
      rmse_train = xgb_metrics_train$rmse,
      mae_train = xgb_metrics_train$mae,
      r2_train = xgb_metrics_train$r2,
      correlation_train = xgb_metrics_train$correlation
    )
  )

  print(results)


  feature_importance <- save_feature_importance(
    models = list(
      lasso = lasso_fit,
      elastic_net = elastic_net_model,
      lightgbm = lgb_model,
      xgboost = xgb_model
    ),
    feature_names = colnames(X_train_scaled),
    model_name = model_name,
    output_dir = model_dir
  )

  return(
    list(
      results = results,
      predictions = predictions,
      train_predictions = train_predictions,
      y_train = y_train,
      y_test = y_test,
      models = list(
        lasso = lasso_fit,
        elastic_net = elastic_net_model,
        lightgbm = lgb_model,
        xgboost = xgb_model
      ),
      feature_importance = feature_importance
    )
  )
}

predict_age <- function(
    model_dir,
    matrix_scaled,
    model_use = c("LASSO", "Elastic Net", "LightGBM", "XGBoost")) {
  model_use <- match.arg(model_use)

  model <- readRDS(model_dir)

  if (model_use %in% c("LASSO", "Elastic Net")) {
    predicted_ages <- predict(model, newx = as.matrix(matrix_scaled))
    if (model_use == "Elastic Net") {
      predicted_ages <- predict(
        model,
        s = model$lambda.min,
        newx = as.matrix(matrix_scaled)
      )
    }
  } else if (model_use == "LightGBM") {
    predicted_ages <- predict(model, newdata = as.matrix(matrix_scaled))
  } else if (model_use == "XGBoost") {
    dtest_xgb <- xgboost::xgb.DMatrix(data = as.matrix(matrix_scaled))
    colnames(dtest_xgb) <- NULL
    predicted_ages <- predict(model, dtest_xgb)
  }

  return(
    data.frame(
      eid = rownames(matrix_scaled),
      predicted_age = as.vector(predicted_ages)
    )
  )
}

predict_ages <- function(
    matrix_scaled,
    model_dir) {
  predictions <- list()

  predictions$lasso <- predict_age(
    file.path(model_dir, "lasso_model.rds"),
    matrix_scaled,
    "LASSO"
  )

  predictions$elastic_net <- predict_age(
    file.path(model_dir, "elastic_net_model.rds"),
    matrix_scaled,
    "Elastic Net"
  )

  predictions$lightgbm <- predict_age(
    file.path(model_dir, "lightgbm_model.rds"),
    matrix_scaled,
    "LightGBM"
  )

  predictions$xgboost <- predict_age(
    file.path(model_dir, "xgboost_model.rds"),
    matrix_scaled,
    "XGBoost"
  )

  pred_df <- data.frame(
    eid = rownames(matrix_scaled),
    lasso_pred = predictions$lasso$predicted_age,
    elastic_net_pred = predictions$elastic_net$predicted_age,
    lightgbm_pred = predictions$lightgbm$predicted_age,
    xgboost_pred = predictions$xgboost$predicted_age
  )

  pred_df$mean_pred <- rowMeans(
    pred_df[, c(
      "lasso_pred", "elastic_net_pred",
      "lightgbm_pred", "xgboost_pred"
    )]
  )

  return(pred_df)
}


save_predictions <- function(
    predictions,
    covariates,
    output_file) {
  results_with_covariates <- merge(
    predictions,
    covariates,
    by = "eid",
    all.x = TRUE
  )

  write.csv(results_with_covariates, output_file, row.names = FALSE)

  return(results_with_covariates)
}

save_feature_importance <- function(
    models,
    feature_names,
    model_name,
    output_dir) {
  model_dir <- file.path(output_dir, model_name)
  if (!dir.exists(model_dir)) {
    dir.create(model_dir, recursive = TRUE)
  }

  lasso_coef <- coef(models$lasso)
  lasso_importance <- data.frame(
    feature = feature_names,
    importance = abs(as.vector(lasso_coef[-1]))
  )
  lasso_importance <- lasso_importance[order(-lasso_importance$importance), ]
  write.csv(
    lasso_importance,
    file.path(model_dir, "lasso_feature_importance.csv"),
    row.names = FALSE
  )

  elastic_net_coef <- coef(models$elastic_net)
  elastic_net_importance <- data.frame(
    feature = feature_names,
    importance = abs(as.vector(elastic_net_coef[-1]))
  )
  elastic_net_importance <- elastic_net_importance[order(-elastic_net_importance$importance), ]
  write.csv(
    elastic_net_importance,
    file.path(model_dir, "elastic_net_feature_importance.csv"),
    row.names = FALSE
  )

  lgb_importance <- lgb.importance(models$lightgbm)
  lgb_importance <- data.frame(
    feature = lgb_importance$Feature,
    importance = lgb_importance$Gain
  )
  lgb_importance <- lgb_importance[order(-lgb_importance$importance), ]
  write.csv(
    lgb_importance,
    file.path(model_dir, "lightgbm_feature_importance.csv"),
    row.names = FALSE
  )

  xgb_importance <- xgb.importance(
    feature_names = feature_names,
    model = models$xgboost
  )
  xgb_importance <- data.frame(
    feature = xgb_importance$Feature,
    importance = xgb_importance$Gain
  )
  xgb_importance <- xgb_importance[order(-xgb_importance$importance), ]
  write.csv(
    xgb_importance,
    file.path(model_dir, "xgboost_feature_importance.csv"),
    row.names = FALSE
  )

  return(
    list(
      lasso = lasso_importance,
      elastic_net = elastic_net_importance,
      lightgbm = lgb_importance,
      xgboost = xgb_importance
    )
  )
}
