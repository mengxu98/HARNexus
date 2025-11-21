create_metric_plots <- function(
    models_results,
    colors_feature_selection = c(
      "HAR" = "#3590bf",
      "HAGR" = "#47a144",
      "BrainAging" = "#de4213",
      "ProtAge" = "#ec651a",
      "SynGO" = "#9947a1"
    ),
    show_value = TRUE,
    add_covariate = TRUE,
    methods_use = c("LASSO", "Elastic Net", "LightGBM", "XGBoost")) {
  all_results <- do.call(
    rbind,
    lapply(
      names(models_results), function(source) {
        results <- models_results[[source]]$results
        results$data_source <- source

        y_test <- models_results[[source]]$y_test
        test_predictions <- models_results[[source]]$predictions

        correlations <- sapply(
          names(test_predictions), function(method) {
            y_pred <- as.vector(test_predictions[[method]])
            calculate_metrics(y_test, y_pred)$correlation
          }
        )

        results$correlation <- correlations[results$method]
        return(results)
      }
    )
  )

  all_results$data_source <- factor(
    all_results$data_source,
    levels = c("HAR", "HAGR", "BrainAging", "ProtAge", "SynGO")
  )

  all_results <- all_results[all_results$method %in% methods_use, ]

  metrics_long <- reshape2::melt(
    all_results,
    id.vars = c("model", "method", "data_source"),
    measure.vars = c("mse", "rmse", "mae", "r2", "correlation"),
    variable.name = "metric",
    value.name = "value"
  )

  metric_plots <- list()

  for (metric_name in c("mse", "rmse", "mae", "r2", "correlation")) {
    metric_data <- metrics_long[metrics_long$metric == metric_name, ]

    p <- ggplot(
      metric_data,
      aes(x = method, y = value, fill = data_source)
    ) +
      geom_bar(
        stat = "identity",
        position = position_dodge(width = 0.9),
        width = 0.8
      ) +
      labs(
        x = "",
        y = ifelse(
          metric_name == "correlation",
          "Correlation",
          ifelse(
            metric_name == "r2",
            expression(R^2),
            toupper(metric_name)
          )
        ),
        fill = "Type"
      ) +
      theme_bw() +
      scale_fill_manual(
        values = colors_feature_selection[levels(metric_data$data_source)],
        breaks = levels(metric_data$data_source)
      ) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))

    if (show_value) {
      p <- p +
        geom_text(
          aes(label = round(value, 2), group = data_source),
          position = position_dodge(width = 0.9),
          vjust = -0.5,
          size = 2.5
        ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
    }

    metric_plots[[metric_name]] <- p
  }

  return(metric_plots)
}

scatter_plot <- function(
    y_train,
    y_test,
    pred_train,
    pred_test,
    method_name,
    source_name,
    color_palette = c(
      "Train set" = "#3366ff",
      "Test set" = "#cc6666"
    ),
    add_covariate = TRUE,
    max_points = 2000,
    use_density = TRUE) {
  cor_train <- cor(y_train, pred_train)
  cor_test <- cor(y_test, pred_test)

  r2_train <- calculate_metrics(y_train, pred_train)$r2
  r2_test <- calculate_metrics(y_test, pred_test)$r2

  train_df <- data.frame(
    Measured = y_train,
    Predicted = pred_train,
    Residuals = pred_train - y_train,
    Set = "Train set"
  )

  test_df <- data.frame(
    Measured = y_test,
    Predicted = pred_test,
    Residuals = pred_test - y_test,
    Set = "Test set"
  )

  combined_df <- rbind(train_df, test_df)

  sample_data <- function(df, group_var, max_per_group) {
    df %>%
      group_by(!!sym(group_var)) %>%
      sample_n(size = min(n(), max_per_group)) %>%
      ungroup()
  }

  if (nrow(combined_df) > max_points) {
    max_per_group <- floor(max_points / length(unique(combined_df$Set)))
    sampled_df <- sample_data(combined_df, "Set", max_per_group)
  } else {
    sampled_df <- combined_df
  }

  p_scatter <- ggplot(combined_df, aes(x = Measured, y = Predicted)) +
    theme_bw()

  if (use_density) {
    p_scatter <- p_scatter +
      stat_density_2d(
        aes(fill = after_stat(density), group = Set),
        geom = "tile",
        contour = FALSE,
        n = 100,
        alpha = 0.7
      ) +
      scale_fill_gradient(low = "white", high = "darkblue", guide = "none")
  }

  p_scatter <- p_scatter +
    geom_point(
      data = sampled_df,
      aes(color = Set),
      alpha = 0.6,
      size = 0.8
    ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      color = "gray45",
      linetype = "solid",
      linewidth = 1
    ) +
    scale_color_manual(values = color_palette) +
    labs(x = "Actual age", y = "Protein predicted age") +
    theme(
      legend.position = "bottom",
      plot.margin = margin(0, 0, 0, 0)
    )

  p_hist_x <- ggplot(combined_df, aes(x = Measured, fill = Set)) +
    geom_density(
      alpha = 0.7,
      position = "stack"
    ) +
    scale_fill_manual(values = color_palette) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0),
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    labs(y = "Density")

  p_hist_y <- ggplot(combined_df, aes(x = Predicted, fill = Set)) +
    geom_density(
      alpha = 0.7,
      position = "stack"
    ) +
    scale_fill_manual(values = color_palette) +
    coord_flip() +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0),
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1)
    ) +
    labs(x = "Density")

  if (nrow(combined_df) > max_points) {
    sampled_df_res <- sample_data(combined_df, "Set", max_per_group)
  } else {
    sampled_df_res <- combined_df
  }

  p_residuals <- ggplot(combined_df, aes(x = Residuals, y = Measured)) +
    theme_bw()

  if (use_density) {
    p_residuals <- p_residuals +
      stat_density_2d(
        aes(fill = after_stat(density), group = Set),
        geom = "tile",
        contour = FALSE,
        n = 100,
        alpha = 0.7
      ) +
      scale_fill_gradient(low = "white", high = "darkred", guide = "none")
  }

  p_residuals <- p_residuals +
    geom_point(
      data = sampled_df_res,
      aes(color = Set),
      alpha = 0.6,
      size = 0.8
    ) +
    geom_vline(
      xintercept = 0,
      color = "gray45",
      linetype = "solid",
      linewidth = 1
    ) +
    scale_color_manual(values = color_palette) +
    theme(
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(angle = 30)
    ) +
    labs(x = "Residuals")

  p_hist_res <- ggplot(combined_df, aes(x = Residuals, fill = Set)) +
    geom_density(
      alpha = 0.7,
      adjust = 1.5,
      position = "stack"
    ) +
    scale_fill_manual(values = color_palette) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    labs(y = "Density")

  design <- "
    AB##
    CDE#
  "

  combined_plot <- p_hist_res + p_hist_x +
    p_residuals + p_scatter + p_hist_y +
    plot_layout(
      widths = c(0.1, 0.3, 0.1),
      heights = c(0.1, 0.4)
    ) +
    plot_layout(design = design) +
    plot_annotation(
      title = paste(
        method_name,
        "\n",
        ifelse(
          add_covariate,
          source_name,
          gsub("_without_covariate", "", source_name)
        )
      ),
      subtitle = paste0(
        " Train correlation: ", round(cor_train, 2),
        " | Test correlation: ", round(cor_test, 2)
      )
    )
  return(combined_plot)
}

create_scatter_plots <- function(
    models_results,
    feature_file,
    base_scatter = TRUE,
    only_test = FALSE,
    color_pattern = c(
      "Train set" = "#3366ff",
      "Test set" = "#cc6666"
    ),
    add_covariate = TRUE,
    max_points = 2000,
    use_density = TRUE) {
  scatter_plots <- list()
  base_scatter_plots <- list()

  for (source in names(models_results)) {
    y_train <- models_results[[source]]$y_train
    y_test <- models_results[[source]]$y_test
    train_predictions <- models_results[[source]]$train_predictions
    test_predictions <- models_results[[source]]$predictions

    for (method_name in names(test_predictions)) {
      y_pred_train <- as.vector(train_predictions[[method_name]])
      y_pred_test <- as.vector(test_predictions[[method_name]])

      if (only_test) {
        test_df <- data.frame(
          Measured = y_test,
          Predicted = y_pred_test,
          Set = "Test set"
        )

        if (nrow(test_df) > max_points) {
          sampled_test_df <- test_df %>%
            sample_n(size = max_points)
        } else {
          sampled_test_df <- test_df
        }

        p_base <- ggplot(test_df, aes(x = Measured, y = Predicted)) +
          theme_bw()

        if (use_density) {
          p_base <- p_base +
            stat_density_2d(
              aes(fill = after_stat(density)),
              geom = "tile",
              contour = FALSE,
              n = 100,
              alpha = 0.7
            ) +
            scale_fill_gradient(low = "white", high = "#cc6666", guide = "none")
        }

        p_base <- p_base +
          geom_point(
            data = sampled_test_df,
            aes(color = Set),
            alpha = 0.7,
            size = 1
          ) +
          geom_abline(
            intercept = 0,
            slope = 1,
            color = "gray45",
            linetype = "solid",
            linewidth = 1
          ) +
          scale_color_manual(values = color_pattern) +
          labs(
            x = "Actual age", y = "Protein predicted age",
            title = paste(
              method_name,
              "\n",
              ifelse(
                add_covariate,
                source_name,
                gsub("_without_covariate", "", source_name)
              )
            ),
            subtitle = paste0(
              "Test correlation: ", round(cor(y_test, y_pred_test), 2)
            )
          ) +
          theme(
            legend.position = "bottom",
            plot.margin = margin(0, 0, 0, 0)
          )

        p_enhanced <- scatter_plot(
          y_train = y_train,
          y_test = y_test,
          pred_train = y_pred_train,
          pred_test = y_pred_test,
          method_name = method_name,
          source_name = source,
          color_palette = color_pattern,
          add_covariate = add_covariate,
          max_points = max_points,
          use_density = use_density
        )

        base_scatter_plots[[paste(source, method_name, sep = "_")]] <- p_base
        scatter_plots[[paste(source, method_name, sep = "_")]] <- p_enhanced
      } else {
        p_enhanced <- scatter_plot(
          y_train = y_train,
          y_test = y_test,
          pred_train = y_pred_train,
          pred_test = y_pred_test,
          method_name = method_name,
          source_name = source,
          color_palette = color_pattern,
          add_covariate = add_covariate,
          max_points = max_points,
          use_density = use_density
        )

        train_df <- data.frame(
          Measured = y_train,
          Predicted = y_pred_train,
          Set = "Train set"
        )

        test_df <- data.frame(
          Measured = y_test,
          Predicted = y_pred_test,
          Set = "Test set"
        )

        combined_df <- rbind(train_df, test_df)

        if (nrow(combined_df) > max_points) {
          max_per_group <- floor(max_points / length(unique(combined_df$Set)))
          sampled_df <- combined_df %>%
            group_by(Set) %>%
            sample_n(size = min(n(), max_per_group)) %>%
            ungroup()
        } else {
          sampled_df <- combined_df
        }

        p_base <- ggplot(combined_df, aes(x = Measured, y = Predicted)) +
          theme_bw()

        if (use_density) {
          p_base <- p_base +
            stat_density_2d(
              aes(fill = after_stat(density), group = Set),
              geom = "tile",
              contour = FALSE,
              n = 100,
              alpha = 0.7
            ) +
            scale_fill_gradient(low = "white", high = "darkblue", guide = "none")
        }

        p_base <- p_base +
          geom_point(
            data = sampled_df,
            aes(color = Set),
            alpha = 0.7,
            size = 1
          ) +
          geom_abline(
            intercept = 0,
            slope = 1,
            color = "gray45",
            linetype = "solid",
            linewidth = 1
          ) +
          scale_color_manual(values = color_pattern) +
          labs(
            x = "Actual age", y = "Protein predicted age",
            title = paste(
              method_name,
              "\n",
              ifelse(
                add_covariate,
                source_name,
                gsub("_without_covariate", "", source_name)
              )
            ),
            subtitle = paste0(
              "Train correlation: ", round(cor(y_train, y_pred_train), 2),
              " | Test correlation: ", round(cor(y_test, y_pred_test), 2)
            )
          ) +
          theme(
            legend.position = "bottom"
          )

        base_scatter_plots[[paste(source, method_name, sep = "_")]] <- p_base
        scatter_plots[[paste(source, method_name, sep = "_")]] <- p_enhanced
      }
    }
  }

  if (base_scatter) {
    return(
      list(
        enhanced_plots = scatter_plots,
        base_plots = base_scatter_plots
      )
    )
  } else {
    return(
      list(
        enhanced_plots = scatter_plots
      )
    )
  }
}


get_pubmed_counts <- function(
    features, pubmed_keywords) {
  pb <- utils::txtProgressBar(min = 0, max = length(features), style = 3)

  counts <- numeric(length(features))
  names(counts) <- features

  if (length(pubmed_keywords) == 1) {
    keywords_query <- pubmed_keywords
  } else {
    keywords_query <- paste(pubmed_keywords, collapse = "+AND+")
  }

  for (i in seq_along(features)) {
    utils::setTxtProgressBar(pb, i)

    feature <- features[i]
    query <- paste0(feature, "+AND+", keywords_query)

    tryCatch(
      {
        url <- paste0(
          "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=",
          query, "&retmode=json"
        )
        response <- httr::GET(url)
        if (httr::status_code(response) == 200) {
          result <- jsonlite::fromJSON(httr::content(response, "text"))

          if (!is.null(result$esearchresult$count)) {
            counts[i] <- as.numeric(result$esearchresult$count)
          }
        }

        Sys.sleep(0.5)
      },
      error = function(e) {
        warning(paste("Error retrieving count for", feature, ":", e$message))
      }
    )
  }

  close(pb)

  result_df <- data.frame(
    Feature = names(counts),
    PubmedCount = as.numeric(counts),
    stringsAsFactors = FALSE
  )

  return(result_df)
}

get_color_gradient <- function(
    count,
    max_count,
    color_palette = c("#62a6b8", "#026077", "gray")) {
  if (is.na(count) || count == 0) {
    return(color_palette[3])
  }

  alpha <- min(1, count / max_count)
  return(
    colorRampPalette(
      color_palette[c(1, 2)]
    )(100)[ceiling(alpha * 99) + 1]
  )
}

create_importance_plot <- function(
    features_df,
    model_name,
    source_name,
    color_palette = c("#62a6b8", "#026077", "gray"),
    add_covariate = TRUE) {
  max_count <- max(features_df$PubmedCount)

  features_df <- features_df %>%
    arrange(desc(Importance)) %>%
    mutate(row_id = row_number())

  features_df$Color <- sapply(
    features_df$PubmedCount, function(count) {
      get_color_gradient(count, max_count, color_palette)
    }
  )

  y_max <- max(features_df$Importance) * 1.2

  p <- ggplot(features_df) +
    geom_bar(
      aes(
        x = reorder(Feature, Importance),
        y = Importance,
        fill = Feature
      ),
      stat = "identity"
    ) +
    scale_fill_manual(
      values = setNames(features_df$Color, features_df$Feature)
    ) +
    geom_line(
      aes(
        x = reorder(Feature, Importance),
        y = PubmedCount / max(PubmedCount) * max(Importance),
        group = 1
      ),
      color = "red", linetype = "dashed"
    ) +
    geom_point(
      aes(
        x = reorder(Feature, Importance),
        y = PubmedCount / max(PubmedCount) * max(Importance)
      ),
      color = "red", size = 1.5
    ) +
    geom_text(
      aes(
        x = reorder(Feature, Importance),
        y = PubmedCount / max(PubmedCount) * max(Importance),
        label = PubmedCount
      ),
      hjust = -0.2, size = 2.5
    ) +
    scale_y_continuous(
      name = "Importance",
      limits = c(0, y_max),
      sec.axis = sec_axis(
        ~ . * max(features_df$PubmedCount) / max(features_df$Importance),
        name = "PubMed counts"
      )
    ) +
    coord_flip() +
    labs(
      title = paste(
        model_name,
        "\n",
        ifelse(
          add_covariate,
          source_name,
          gsub("_without_covariate", "", source_name)
        )
      )
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title.y = element_blank(),
      plot.margin = margin(5, 20, 5, 5)
    )

  return(p)
}

analyze_feature_importance <- function(
    models_results,
    top_n = 10,
    color_palette = c("#62a6b8", "#026077", "gray"),
    pubmed_counts_file = NULL,
    pubmed_keyword = NULL,
    add_covariate = TRUE) {
  all_important_features <- list()

  for (source in names(models_results)) {
    lasso_model <- models_results[[source]]$models$lasso
    lasso_coef <- as.matrix(coef(lasso_model))
    lasso_coef <- lasso_coef[-1, , drop = FALSE]

    lasso_coef <- lasso_coef[rownames(lasso_coef) != "gender", , drop = FALSE]
    lasso_features <- rownames(lasso_coef)[lasso_coef != 0]

    elastic_net_model <- models_results[[source]]$models$elastic_net
    elastic_net_coef <- as.matrix(
      coef(elastic_net_model, s = elastic_net_model$lambda.min)
    )
    elastic_net_coef <- elastic_net_coef[-1, , drop = FALSE]
    elastic_net_coef <- elastic_net_coef[rownames(elastic_net_coef) != "gender", , drop = FALSE]
    elastic_net_features <- rownames(elastic_net_coef)[elastic_net_coef != 0]

    lgb_model <- models_results[[source]]$models$lightgbm
    lgb_importance <- lightgbm::lgb.importance(lgb_model)
    lgb_importance <- lgb_importance[lgb_importance$Feature != "gender", ]
    lgb_features <- lgb_importance$Feature

    xgb_model <- models_results[[source]]$models$xgboost
    xgb_importance <- xgboost::xgb.importance(model = xgb_model)
    xgb_importance <- xgb_importance[xgb_importance$Feature != "gender", ]
    xgb_features <- xgb_importance$Feature

    all_important_features[[source]] <- list(
      lasso = lasso_features,
      elastic_net = elastic_net_features,
      lightgbm = lgb_features,
      xgboost = xgb_features
    )
  }

  unique_features <- unique(
    unlist(
      lapply(
        all_important_features, function(x) {
          unlist(x)
        }
      )
    )
  )

  log_message(sprintf("Total unique features: %d\n", length(unique_features)))

  if (!is.null(pubmed_counts_file)) {
    if (file.exists(pubmed_counts_file)) {
      pubmed_data <- read.csv(pubmed_counts_file)
      pubmed_counts <- setNames(pubmed_data$PubmedCount, pubmed_data$Feature)
      log_message("Loaded PubMed counts from file.\n")
    } else {
      log_message("File does not exist. Getting PubMed counts...")
      pubmed_counts <- get_pubmed_counts(
        features = unique_features,
        pubmed_keywords = pubmed_keyword
      )
      write.csv(
        pubmed_counts,
        pubmed_counts_file,
        row.names = FALSE
      )
    }
  } else {
    log_message("Getting PubMed counts...")
    pubmed_counts <- get_pubmed_counts(
      features = unique_features,
      pubmed_keywords = pubmed_keyword
    )
    write.csv(
      pubmed_counts,
      "feature_pubmed_counts.csv",
      row.names = FALSE
    )
  }

  importance_plots <- list()

  for (source in names(models_results)) {
    lasso_model <- models_results[[source]]$models$lasso
    lasso_coef <- as.matrix(coef(lasso_model))
    lasso_coef <- lasso_coef[-1, , drop = FALSE]

    lasso_coef <- lasso_coef[rownames(lasso_coef) != "gender", , drop = FALSE]

    non_zero_coef <- lasso_coef[lasso_coef != 0, , drop = FALSE]
    if (length(non_zero_coef) > 0) {
      important_features <- data.frame(
        Feature = rownames(lasso_coef)[lasso_coef != 0],
        Importance = abs(non_zero_coef)[, 1],
        Model = "LASSO",
        Dataset = source
      )

      important_features$PubmedCount <- sapply(
        important_features$Feature, function(f) {
          ifelse(f %in% names(pubmed_counts), pubmed_counts[f], 0)
        }
      )

      if (nrow(important_features) > top_n) {
        important_features$Importance <- as.numeric(important_features$Importance)
        important_features$Importance[is.na(important_features$Importance)] <- 0
        important_features <- important_features[order(
          important_features$Importance,
          decreasing = TRUE
        ), ][1:top_n, ]
      }

      p <- create_importance_plot(
        features_df = important_features,
        model_name = "LASSO",
        source_name = source,
        color_palette = color_palette,
        add_covariate = add_covariate
      )
      importance_plots[[paste(source, "LASSO", sep = "_")]] <- p
    }

    elastic_net_model <- models_results[[source]]$models$elastic_net
    elastic_net_coef <- as.matrix(
      coef(elastic_net_model,
        s = elastic_net_model$lambda.min
      )
    )
    elastic_net_coef <- elastic_net_coef[-1, , drop = FALSE]

    elastic_net_coef <- elastic_net_coef[rownames(elastic_net_coef) != "gender", , drop = FALSE]

    non_zero_coef <- elastic_net_coef[elastic_net_coef != 0, , drop = FALSE]
    if (length(non_zero_coef) > 0) {
      important_features <- data.frame(
        Feature = rownames(elastic_net_coef)[elastic_net_coef != 0],
        Importance = abs(non_zero_coef)[, 1],
        Model = "Elastic Net",
        Dataset = source
      )

      important_features$PubmedCount <- sapply(
        important_features$Feature, function(f) {
          ifelse(f %in% names(pubmed_counts), pubmed_counts[f], 0)
        }
      )

      if (nrow(important_features) > top_n) {
        important_features$Importance <- as.numeric(important_features$Importance)
        important_features$Importance[is.na(important_features$Importance)] <- 0
        important_features <- important_features[order(
          important_features$Importance,
          decreasing = TRUE
        ), ][1:top_n, ]
      }

      p <- create_importance_plot(
        features_df = important_features,
        model_name = "Elastic Net",
        source_name = source,
        color_palette = color_palette,
        add_covariate = add_covariate
      )
      importance_plots[[paste(source, "Elastic Net", sep = "_")]] <- p
    }

    lgb_model <- models_results[[source]]$models$lightgbm
    lgb_importance <- lgb.importance(lgb_model)

    lgb_importance <- lgb_importance[lgb_importance$Feature != "gender", ]

    if (nrow(lgb_importance) > 0) {
      if (nrow(lgb_importance) > top_n) {
        lgb_importance <- lgb_importance[1:min(top_n, nrow(lgb_importance)), ]
      }

      lgb_importance_df <- data.frame(
        Feature = lgb_importance$Feature,
        Importance = lgb_importance$Gain,
        Model = "LightGBM",
        Dataset = source
      )

      lgb_importance_df$PubmedCount <- sapply(
        lgb_importance_df$Feature, function(f) {
          ifelse(f %in% names(pubmed_counts), pubmed_counts[f], 0)
        }
      )

      p <- create_importance_plot(
        features_df = lgb_importance_df,
        model_name = "LightGBM",
        source_name = source,
        color_palette = color_palette,
        add_covariate = add_covariate
      )
      importance_plots[[paste(source, "LightGBM", sep = "_")]] <- p
    }

    xgb_model <- models_results[[source]]$models$xgboost
    xgb_importance <- xgboost::xgb.importance(model = xgb_model)

    xgb_importance <- xgb_importance[xgb_importance$Feature != "gender", ]

    if (nrow(xgb_importance) > 0) {
      if (nrow(xgb_importance) > top_n) {
        xgb_importance <- xgb_importance[1:min(top_n, nrow(xgb_importance)), ]
      }

      xgb_importance_df <- data.frame(
        Feature = xgb_importance$Feature,
        Importance = xgb_importance$Gain,
        Model = "XGBoost",
        Dataset = source
      )

      xgb_importance_df$PubmedCount <- sapply(
        xgb_importance_df$Feature, function(f) {
          ifelse(f %in% names(pubmed_counts), pubmed_counts[f], 0)
        }
      )

      p <- create_importance_plot(
        features_df = xgb_importance_df,
        model_name = "XGBoost",
        source_name = source,
        color_palette = color_palette,
        add_covariate = add_covariate
      )
      importance_plots[[paste(source, "XGBoost", sep = "_")]] <- p
    }
  }

  return(importance_plots)
}

create_performance_table <- function(
    models_results,
    path = "performance_table.docx") {
  results_table <- do.call(
    rbind,
    lapply(
      names(models_results), function(source) {
        models_results[[source]]$results$data_source <- source
        return(models_results[[source]]$results)
      }
    )
  )

  results_table$model <- gsub(
    "HAR_without_covariate",
    "HAR aging model",
    results_table$model
  )
  results_table$model <- gsub(
    "HAR_with_covariate",
    "HAR aging model (+gender)",
    results_table$model
  )
  results_table$model <- gsub(
    "HAR",
    "HAR aging model",
    results_table$model
  )
  results_table$model <- gsub(
    "HAGR",
    "HAGR aging model",
    results_table$model
  )
  results_table$model <- gsub(
    "HAGR_without_covariate",
    "HAGR aging model",
    results_table$model
  )
  results_table$model <- gsub(
    "HAGR_with_covariate",
    "HAGR aging model (+gender)",
    results_table$model
  )
  results_table$model <- gsub(
    "ProtAge",
    "ProtAge model",
    results_table$model
  )
  results_table$model <- gsub(
    "ProtAge_without_covariate",
    "ProtAge model",
    results_table$model
  )

  table_data <- results_table[, c(
    "model", "method", "mse", "rmse", "mae", "r2", "correlation"
  )]
  table_data$mse <- round(table_data$mse, 2)
  table_data$rmse <- round(table_data$rmse, 2)
  table_data$mae <- round(table_data$mae, 2)
  table_data$r2 <- round(table_data$r2, 2)

  ft <- flextable(table_data)

  ft <- set_caption(
    ft,
    caption = "Table *. Performance comparison of different models"
  )

  ft <- set_header_labels(ft,
    model = "Model",
    method = "Method",
    mse = "MSE",
    rmse = "RMSE",
    mae = "MAE",
    r2 = "R²",
    correlation = "Correlation"
  )

  ft <- bold(ft, part = "header")

  best_model_mse <- table_data$method[which.min(table_data$mse)]
  best_model_rmse <- table_data$method[which.min(table_data$rmse)]
  best_model_mae <- table_data$method[which.min(table_data$mae)]
  best_model_r2 <- table_data$method[which.max(table_data$r2)]
  best_model_correlation <- table_data$method[which.max(table_data$correlation)]

  ft <- bold(ft,
    i = which(table_data$method == best_model_mse),
    j = c("mse")
  )
  ft <- bold(ft,
    i = which(table_data$method == best_model_rmse),
    j = c("rmse")
  )
  ft <- bold(ft,
    i = which(table_data$method == best_model_mae),
    j = c("mae")
  )
  ft <- bold(ft,
    i = which(table_data$method == best_model_r2),
    j = c("r2")
  )
  ft <- bold(ft,
    i = which(table_data$method == best_model_correlation),
    j = c("correlation")
  )


  ft <- border(
    ft,
    part = "header",
    border.top = fp_border(width = 2)
  )
  ft <- border(
    ft,
    part = "header",
    border.bottom = fp_border(width = 1)
  )
  ft <- border(
    ft,
    part = "body",
    border.bottom = fp_border(width = 2),
    i = nrow(table_data)
  )

  ft <- theme_booktabs(ft)
  ft <- autofit(ft)

  doc <- read_docx()
  doc <- body_add_flextable(doc, value = ft)
  print(doc, target = path)

  return(ft)
}


create_venn_diagram <- function(
    feature_files,
    colors_feature_selection = c(
      "HAR" = "#3590bf",
      "HAGR" = "#47a144",
      "BrainAging" = "#de4213",
      "ProtAge" = "#ec651a",
      "SynGO" = "#9947a1"
    )) {
  genes_set_names <- names(feature_files)
  genes_set <- lapply(genes_set_names, function(genes_set_name) {
    read.csv(feature_files[[genes_set_name]])[, 1]
  })
  names(genes_set) <- genes_set_names

  venn_plot <- ggVennDiagram::ggVennDiagram(
    genes_set,
    category.names = names(genes_set),
    label_alpha = 0,
    label_color = "white",
    label_size = 5,
    set_color = colors_feature_selection[genes_set_names],
    edge_lty = "solid",
    edge_size = 1.5
  )

  return(venn_plot)
}

create_citation_boxplot <- function(
    feature_files,
    pubmed_counts_file,
    colors_feature_selection = c(
      "HAR" = "#3590bf",
      "HAGR" = "#47a144",
      "BrainAging" = "#de4213",
      "ProtAge" = "#ec651a",
      "SynGO" = "#9947a1"
    ),
    add_signif = FALSE) {
  pubmed_data <- read.csv(pubmed_counts_file)

  genes_set_names <- names(feature_files)
  genes_set <- lapply(genes_set_names, function(genes_set_name) {
    read.csv(feature_files[[genes_set_name]])[, 1]
  })
  names(genes_set) <- genes_set_names

  har_genes <- genes_set[["HAR"]]
  hagr_genes <- genes_set[["HAGR"]]

  plot_data <- data.frame(
    Gene = c(har_genes, hagr_genes),
    Group = factor(
      c(
        rep("HAR", length(har_genes)),
        rep("HAGR", length(hagr_genes))
      ),
      levels = c("HAR", "HAGR")
    )
  )

  plot_data$counts <- sapply(
    plot_data$Gene, function(g) {
      idx <- which(pubmed_data$Feature == g)
      if (length(idx) > 0) pubmed_data$PubmedCount[idx] else 0
    }
  )

  plot_data$counts <- log1p(plot_data$counts)

  citation_plot <- ggplot(
    plot_data,
    aes(x = Group, y = counts, fill = Group)
  ) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) +
    scale_fill_manual(
      values = colors_feature_selection[unique(plot_data$Group)]
    ) +
    labs(y = "log10(PubMed counts + 1)", title = "PubMed counts") +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = 20)
    )

  if (add_signif) {
    groups <- levels(plot_data$Group)
    test_results <- data.frame()
    for (i in 1:(length(groups) - 1)) {
      for (j in (i + 1):length(groups)) {
        group1 <- groups[i]
        group2 <- groups[j]
        test_data <- plot_data[plot_data$Group %in% c(group1, group2), ]
        test_result <- suppressWarnings(
          wilcox.test(
            counts ~ Group,
            data = test_data
          )
        )

        test_results <- rbind(test_results, data.frame(
          Comparison = paste(group1, "vs", group2),
          Group1 = group1,
          Group2 = group2,
          W_statistic = test_result$statistic,
          P_value = test_result$p.value,
          Significant = ifelse(test_result$p.value < 0.05, "Yes", "No")
        ))
      }
    }

    sig_comparisons <- sum(test_results$P_value < 0.05)

    y_max <- max(plot_data$counts)
    y_buffer <- y_max * 0.5

    y_step <- y_buffer / (sig_comparisons + 1)
    y_pos <- y_max + y_step

    print(test_results)

    significant_idx <- which(test_results$P_value < 0.05)
    for (i in seq_along(significant_idx)) {
      idx <- significant_idx[i]

      stars <- ifelse(test_results$P_value[idx] < 0.001, "***",
        ifelse(test_results$P_value[idx] < 0.01, "**",
          ifelse(test_results$P_value[idx] < 0.05, "*", "")
        )
      )

      group1_pos <- which(groups == test_results$Group1[idx])
      group2_pos <- which(groups == test_results$Group2[idx])

      citation_plot <- citation_plot +
        geom_segment(
          aes(
            x = group1_pos,
            xend = group2_pos,
            y = y_max + i * y_step,
            yend = y_max + i * y_step
          ),
          color = "black",
          linewidth = 0.7
        ) +
        geom_text(
          aes(
            x = (group1_pos + group2_pos) / 2,
            y = y_max + i * y_step + y_step * 0.3
          ),
          label = stars,
          size = 6
        )
    }

    citation_plot <- citation_plot +
      scale_y_continuous(
        limits = c(0, y_max + (sig_comparisons + 1) * y_step),
        expand = expansion(mult = c(0, 0.1))
      )
  }

  return(
    list(
      plot = citation_plot,
      test_results = if (add_signif) test_results else NULL
    )
  )
}
