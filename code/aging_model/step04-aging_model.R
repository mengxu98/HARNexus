source("code/functions/aging_data.R")
source("code/functions/aging_model.R")
source("code/functions/aging_plots.R")
source("code/functions/packages.R")
source("code/functions/utils.R")

option_list <- list(
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
    default = "results/aging_model_60_100/models/",
    help = "Base directory for results [default %default]",
    metavar = "character"
  ),
  make_option(c("-e", "--re_train_models"),
    type = "logical",
    default = FALSE,
    help = "Whether to re-train models [default %default]",
    metavar = "logical"
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
  make_option(c("-f", "--feature_files"),
    type = "character",
    default = "HAR:results/aging_model_20_100/features_selection/har/features_optimized.csv,HAGR:results/aging_model_20_100/features_selection/hagr/features_optimized.csv",
    help = "Comma-separated list of feature files in format 'name:path'. [default %default]",
    metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

set.seed(opt$seed)

res_dir <- check_dir(opt$output_dir)
fig_dir <- check_dir(paste0(res_dir, "/plots"))
fig_scatter_dir1 <- check_dir(paste0(fig_dir, "/scatter plots"))
fig_scatter_dir2 <- check_dir(paste0(fig_dir, "/scatter plots_base"))

feature_files <- list()
if (!is.null(opt$feature_files)) {
  file_pairs <- strsplit(opt$feature_files, ",")[[1]]
  for (pair in file_pairs) {
    parts <- strsplit(pair, ":")[[1]]
    if (length(parts) == 2) {
      feature_files[[parts[1]]] <- parts[2]
    }
  }
}

feature_files <- list(
  HAR = "results/aging_model_60_100/features_selection/har/features_optimized.csv",
  HAGR = "results/aging_model_60_100/features_selection/hagr/features_optimized.csv",
  BrainAging = "results/aging_model_60_100/features_selection/brain/features_optimized.csv"
)

data_file <- opt$data_file
response_file <- opt$response_file
results_file <- paste0(res_dir, "/models_results_", opt$train_ratio, ".rds")

if (opt$re_train_models) {
  log_message("Removing existing models...")
  unlink(results_file)
}

if (!file.exists(results_file)) {
  models_results <- list()
  for (source_name in names(feature_files)) {
    log_message(sprintf("Using %s Dataset", source_name))

    data <- prepare_data_new(
      data_file = data_file,
      feature_file = feature_files[[source_name]],
      response_file = response_file,
      preprocessing = FALSE,
      precentage = opt$train_ratio,
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

    models_results[[source_name]] <- train_evaluate_models(
      data_sets$X_train_scaled,
      data_sets$X_test_scaled,
      data_sets$y_train,
      data_sets$y_test,
      model_name = source_name,
      model_dir = paste0(res_dir, "/models/")
    )
  }
  saveRDS(
    models_results,
    file = results_file
  )
} else {
  models_results <- readRDS(results_file)
}

log_message("Creating Performance Table")
performance_table <- create_performance_table(
  models_results,
  path = paste0(res_dir, "/performance_table.docx")
)

log_message("Creating Result metric_plots_list")
metric_plots_list <- create_metric_plots(
  models_results = models_results,
  add_covariate = FALSE,
  methods_use = c("LASSO", "Elastic Net", "LightGBM")
)


log_message("Analyzing Feature Importance")
importance_plots <- analyze_feature_importance(
  models_results,
  top_n = 10,
  pubmed_counts_file = paste0(
    res_dir,
    "/feature_pubmed_counts_brain_aging.csv"
  ),
  pubmed_keyword = c("Brain", "Aging"),
  color_palette = c("#709dff", "#0537a3", "gray"),
  add_covariate = FALSE
)

fig_s2a <- metric_plots_list[[1]] +
  metric_plots_list[[2]] +
  metric_plots_list[[3]] +
  metric_plots_list[[4]] +
  metric_plots_list[[5]] +
  patchwork::plot_layout(
    ncol = 5,
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom"
  )
ggsave(
  paste0(fig_dir, "/Fig.6b.pdf"),
  fig_s2a,
  width = 10,
  height = 3.5
)


fig_s2b <-
  importance_plots[[2]] +
  importance_plots[[6]] +
  importance_plots[[10]] +
  importance_plots[[1]] +
  importance_plots[[5]] +
  importance_plots[[9]] +
  patchwork::plot_layout(
    ncol = 3,
    guides = "collect"
  )

ggsave(
  paste0(fig_dir, "/Fig.S9.pdf"),
  fig_s2b,
  width = 8,
  height = 6
)

log_message("Creating Gene Sets Venn Diagram")
venn_plot <- create_venn_diagram(feature_files)
ggsave(
  paste0(fig_dir, "/fig.6a.pdf"),
  venn_plot,
  width = 3.5,
  height = 3.5
)

citation_comparison <- create_citation_boxplot(
  feature_files,
  pubmed_counts_file = paste0(
    res_dir, "/feature_pubmed_counts_brain_aging.csv"
  ),
  add_signif = TRUE
)
citation_test <- citation_comparison$test_results
write.csv(
  citation_test,
  file = paste0(res_dir, "/citation_test.csv"),
  row.names = FALSE
)
citation_plot <- citation_comparison$plot
citation_plot <- citation_plot +
  labs(title = "PubMed counts", y = "log10(counts)") +
  theme(
    legend.position = "right"
  )
ggsave(
  paste0(fig_dir, "/fig.6c.pdf"),
  citation_plot,
  width = 2.5,
  height = 2.5
)

fig_6_de <- importance_plots[[3]] +
  importance_plots[[7]] +
  importance_plots[[11]] +
  plot_layout(ncol = 3)
ggsave(
  paste0(fig_dir, "/fig.6de.pdf"),
  fig_6_de,
  width = 10,
  height = 3.5
)


log_message("Saving Scatter Plots")
scatter_plot_list <- create_scatter_plots(
  models_results,
  feature_files,
  add_covariate = FALSE,
  max_points = 1000,
  use_density = TRUE
)
scatter_enhanced_plot_list <- scatter_plot_list$enhanced_plots
for (plot_name in names(scatter_enhanced_plot_list)) {
  log_message(sprintf("Saving %s...", plot_name))

  ggsave(
    paste0(
      fig_scatter_dir1, "/", plot_name, ".pdf"
    ),
    scatter_enhanced_plot_list[[plot_name]],
    width = 6,
    height = 5.5
  )
}

scatter_base_plot_list <- scatter_plot_list$base_plots
for (plot_name in names(scatter_base_plot_list)) {
  log_message(sprintf("Saving %s...", plot_name))
  ggsave(
    paste0(
      fig_scatter_dir2, "/", plot_name, "_base.pdf"
    ),
    scatter_base_plot_list[[plot_name]],
    width = 3.5,
    height = 4.5
  )
}
