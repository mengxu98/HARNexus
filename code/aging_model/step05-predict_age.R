source("code/functions/aging_model.R")

if (!file.exists("data/ukb/olink_data_imputed_mean_drop30_har_genes.csv")) {
  expression_matrix <- read.csv(
    "data/ukb/olink_data_imputed_mean_drop30.csv",
    row.names = 1
  )

  features_selection <- read.csv(
    "results/aging_model/features_selection/har/features_optimized.csv"
  )[, 1]

  expression_matrix <- expression_matrix[, features_selection]
  write.csv(
    expression_matrix,
    "data/ukb/olink_data_imputed_mean_drop30_har_genes.csv",
    row.names = TRUE
  )
} else {
  expression_matrix <- read.csv(
    "data/ukb/olink_data_imputed_mean_drop30_har_genes.csv",
    row.names = 1
  )
}

covariates <- read.csv(
  "data/ukb/granular_age_april_07_2025.csv",
  row.names = 1
)


covariates$eid <- rownames(covariates)
common_eids <- intersect(rownames(expression_matrix), covariates$eid)
expression_matrix <- expression_matrix[common_eids, ]
covariates <- covariates[common_eids, ]

predictions <- predict_ages(
  matrix_scaled = scale(expression_matrix),
  model_dir = "results/aging_model_20_100/models/models/HAR/"
)

final_results <- save_predictions(
  predictions = predictions,
  covariates = covariates,
  output_file = "results/aging_model_20_100/models/predicted_ages_HAR.csv"
)

predicted_ages <- read.csv(
  paste0("results/aging_model_20_100/models/predicted_ages_HAR.csv")
)
predicted_ages1 <- read.csv(
  "predicted_ages_without_covariate_HAR.csv"
)


cor(predicted_ages$elastic_net_pred, predicted_ages$age_granular)
cor(predicted_ages$lightgbm_pred, predicted_ages$age_granular)
cor(predicted_ages$xgboost_pred, predicted_ages$age_granular)
cor(predicted_ages$lasso_pred, predicted_ages$age_granular)
cor(predicted_ages1$lightgbm_pred, predicted_ages1$precise_age_at_recruitment)
cor(predicted_ages1$xgboost_pred, predicted_ages1$precise_age_at_recruitment)
cor(predicted_ages1$lasso_pred, predicted_ages1$precise_age_at_recruitment)
cor(predicted_ages1$elastic_net_pred, predicted_ages1$precise_age_at_recruitment)


cor_data <- final_results[, c(
  "lasso_pred", "elastic_net_pred",
  "lightgbm_pred", "xgboost_pred"
)]
colnames(cor_data) <- c("LASSO", "Elastic Net", "LightGBM", "XGBoost")
model_correlations <- cor(cor_data)
log_message("\nCorrelation between models:\n")
print(model_correlations)
pdf(
  "results/aging_model/models/plots/model_correlations.pdf",
  width = 4.2,
  height = 3.2
)
ComplexHeatmap::Heatmap(
  model_correlations,
  name = "Correlation",
  col = colorRamp2(c(0.8, 1), c("#ffffff", "#b61313")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      sprintf("%.2f", model_correlations[i, j]),
      x, y,
      gp = gpar(fontsize = 12)
    )
  },
  heatmap_width = unit(1, "npc"),
  width = 3,
  heatmap_height = unit(1, "npc"),
  height = 3,
  rect_gp = gpar(col = "white", lwd = 1.5),
  border = TRUE,
  column_names_rot = 45
)
dev.off()
