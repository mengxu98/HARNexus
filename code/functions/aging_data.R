source("code/functions/packages.R")

process_protein_data <- function(
    input_path,
    output_path,
    imputation = TRUE,
    zero_imputation = FALSE,
    remove_missing_values = FALSE) {
  log_message("Processing protein data...\n")

  data <- read.csv(input_path, row.names = 1)
  rownames_data <- rownames(data)

  if (zero_imputation) {
    log_message("Processing zero imputation...")
    print(dim(data))
    data <- as.matrix(data)
    data[is.na(data)] <- 0
    write.csv(data, output_path)
    print(dim(data))
    log_message("Zero imputation completed!\n")
    return(data)
  }

  if (imputation) {
    log_message("Imputing missing values...")
    missing_prop <- colMeans(is.na(data))

    valid_columns <- names(missing_prop[missing_prop < 0.3])
    data <- data[, valid_columns]

    log_message("Total number of proteins:", ncol(data), "\n")

    log_message("Starting imputation...\n")
    imputed_data <- mice(
      data,
      m = 3,
      maxit = 3,
      method = "pmm",
      seed = 42,
      print = FALSE
    )

    complete_data <- complete(imputed_data, 1)
    rownames(complete_data) <- rownames_data

    write.csv(complete_data, output_path)
    log_message("Imputation completed!\n")

    return(complete_data)
  } else {
    log_message("Removing missing values...")
    data <- na.omit(data)
    write.csv(data, output_path)
    log_message("Missing values removed!\n")
    return(data)
  }
}

process_participant_data <- function(
    input_path,
    output_path) {
  log_message("Processing participant data...\n")

  data <- read.csv(input_path)
  data$Month.of.birth <- match(data$p52, month.name)
  data$Year.of.birth <- as.numeric(data$p34)
  data$birth_date <- as.Date(
    paste(
      data$p34,
      sprintf("%02d", data$Month.of.birth),
      "01",
      sep = "-"
    )
  )

  data$recruitment_date <- data$birth_date + years(data$p21022)

  data$age_granular <- as.numeric(
    difftime(
      data$recruitment_date,
      data$birth_date,
      units = "days"
    )
  ) / 365.25

  write.csv(data, output_path, row.names = FALSE)
  log_message("Participant data processing completed!\n")

  return(data)
}

preprocess_data <- function(
    X,
    y,
    gender_numeric,
    remove_outliers = TRUE,
    outliers_threshold = 2,
    quantile_threshold_lower = 0.10,
    quantile_threshold_upper = 0.90) {
  if (remove_outliers) {
    outliers <- logical(nrow(X))
    for (col in colnames(X)) {
      Q1 <- quantile(X[[col]], quantile_threshold_lower)
      Q3 <- quantile(X[[col]], quantile_threshold_upper)
      IQR <- Q3 - Q1
      lower_bound <- Q1 - outliers_threshold * IQR
      upper_bound <- Q3 + outliers_threshold * IQR
      outliers <- outliers | (X[[col]] < lower_bound | X[[col]] > upper_bound)
    }

    X <- X[!outliers, ]
    y <- y[!outliers]
    gender_numeric <- gender_numeric[!outliers]

    log_message(sprintf("Removed %d samples with extreme values\n", sum(outliers)))
  }

  return(
    list(
      X = X,
      y = y,
      gender_numeric = gender_numeric
    )
  )
}

prepare_data <- function(
    feature_file,
    response_file,
    samples = NULL,
    preprocessing = FALSE,
    precentage = 0.7,
    age_at = 40,
    use_pca = FALSE) {
  X <- read.csv(feature_file, row.names = 1)
  colnames(X) <- toupper(colnames(X))

  y <- read.csv(response_file)
  y <- y[y$age_granular >= age_at, ]

  common_eids <- intersect(rownames(X), y$eid)
  if (!is.null(samples)) {
    common_eids <- intersect(common_eids, samples)
  }
  X <- X[common_eids, ]
  y <- y[y$eid %in% common_eids, ]

  X <- X[order(rownames(X)), ]
  y <- y[order(y$eid), ]
  response_variable <- y$age_granular
  stopifnot(all(rownames(X) == y$eid))

  gender_numeric <- ifelse(y$p31 == "Male", 1, 0)

  if (preprocessing) {
    log_message("Performing enhanced preprocessing...\n")
    preprocessed_data <- preprocess_data(
      X,
      response_variable,
      gender_numeric,
      remove_outliers = TRUE
    )
    X <- preprocessed_data$X
    response_variable <- preprocessed_data$y
    gender_numeric <- preprocessed_data$gender_numeric
  }
  if (use_pca) {
    pca_model <- prcomp(X, scale. = TRUE)
    X <- predict(pca_model, X)
  }

  set.seed(123)
  train_index <- createDataPartition(
    response_variable,
    p = precentage,
    list = FALSE
  )

  log_message(
    sprintf(
      "Final feature matrix dimensions: %d x %d\n",
      nrow(X), ncol(X)
    )
  )

  return(
    list(
      X = X,
      response_variable = response_variable,
      gender_numeric = gender_numeric,
      train_index = train_index
    )
  )
}

prepare_data_new <- function(
    data_file,
    feature_file = NULL,
    response_file,
    samples = NULL,
    preprocessing = FALSE,
    precentage = 0.7,
    age_at = 40,
    use_pca = FALSE) {
  X <- read.csv(data_file, row.names = 1)
  colnames(X) <- toupper(colnames(X))
  log_message(sprintf("Before filtering: %d samples...", nrow(X)))
  log_message(sprintf("Before filtering: %d features...", ncol(X)))

  if (!is.null(feature_file)) {
    features <- read.csv(feature_file)[, 1]
    features <- toupper(features)
    features <- intersect(colnames(X), features)
    log_message(sprintf("After filtering: %d features...", length(features)))
    X <- X[, features]
  }

  y <- read.csv(response_file)
  y <- y[y$age_granular >= age_at, ]

  common_eids <- intersect(rownames(X), y$eid)
  if (!is.null(samples)) {
    common_eids <- intersect(common_eids, samples)
  }
  X <- X[common_eids, ]
  y <- y[y$eid %in% common_eids, ]

  X <- X[order(rownames(X)), ]
  y <- y[order(y$eid), ]
  response_variable <- y$age_granular
  stopifnot(all(rownames(X) == y$eid))

  gender_numeric <- ifelse(y$p31 == "Male", 1, 0)

  if (preprocessing) {
    log_message("Performing enhanced preprocessing...\n")
    preprocessed_data <- preprocess_data(
      X,
      response_variable,
      gender_numeric,
      remove_outliers = TRUE
    )
    X <- preprocessed_data$X
    response_variable <- preprocessed_data$y
    gender_numeric <- preprocessed_data$gender_numeric
  }
  if (use_pca) {
    pca_model <- prcomp(X, scale. = TRUE)
    X <- predict(pca_model, X)
  }

  set.seed(123)
  train_index <- createDataPartition(
    response_variable,
    p = precentage,
    list = FALSE
  )

  log_message(
    sprintf(
      "Final feature matrix dimensions: %d x %d\n",
      nrow(X), ncol(X)
    )
  )

  return(
    list(
      X = X,
      response_variable = response_variable,
      gender_numeric = gender_numeric,
      train_index = train_index
    )
  )
}

create_train_test_sets <- function(
    X,
    response_variable,
    gender_numeric,
    train_index,
    add_covariate = FALSE,
    scale = TRUE) {
  X_train <- X[train_index, ]
  X_test <- X[-train_index, ]
  y_train <- response_variable[train_index]
  y_test <- response_variable[-train_index]

  if (add_covariate) {
    gender_train <- gender_numeric[train_index]
    gender_test <- gender_numeric[-train_index]

    X_train$gender <- gender_train
    X_test$gender <- gender_test

    if (scale) {
      X_train_scaled <- scale(X_train)
      X_test_scaled <- scale(X_test,
        center = attr(X_train_scaled, "scaled:center"),
        scale = attr(X_train_scaled, "scaled:scale")
      )
    } else {
      X_train_scaled <- X_train
      X_test_scaled <- X_test
    }
  } else {
    if (scale) {
      X_train_scaled <- scale(
        X_train
      )
      X_test_scaled <- scale(
        X_test,
        center = attr(X_train_scaled, "scaled:center"),
        scale = attr(X_train_scaled, "scaled:scale")
      )
    } else {
      X_train_scaled <- X_train
      X_test_scaled <- X_test
    }
  }

  return(
    list(
      X_train_scaled = X_train_scaled,
      X_test_scaled = X_test_scaled,
      y_train = y_train,
      y_test = y_test
    )
  )
}
