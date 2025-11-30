source("code/functions/prepare_env.R")

data_dir <- "../../data/BrainData/integration/"
fig_dir <- check_dir("results/figures/lisi")

if (!file.exists(paste0(data_dir, "/lisi_results.rds"))) {
  lisi_data <- readRDS(paste0(data_dir, "/lisi_data.rds"))
  before_lisi <- compute_lisi(
    X = lisi_data[[3]],
    meta_data = lisi_data[[5]],
    label_colnames = "Dataset"
  )
  after_lisi <- compute_lisi(
    X = lisi_data[[4]],
    meta_data = lisi_data[[5]],
    label_colnames = "Dataset"
  )
  lisi_results <- cbind(
    before_lisi,
    after_lisi
  )
  names(lisi_results) <- c("Raw", "RPCA")

  saveRDS(lisi_results, file.path(data_dir, "lisi_results.rds"))
} else {
  lisi_results <- readRDS(file.path(data_dir, "lisi_results.rds"))
}

lisi_long <- tidyr::gather(
  lisi_results,
  key = "Method",
  value = "LISI"
)
lisi_long$Method <- factor(
  lisi_long$Method,
  levels = c("Raw", "RPCA")
)

p <- ggplot(lisi_long, aes(x = Method, y = LISI)) +
  geom_boxplot(
    aes(fill = Method),
    # alpha = 0.7,
    width = 0.5,
    outlier.shape = NA
  ) +
  scale_fill_manual(
    values = c("Raw" = "#E41A1C", "RPCA" = "#377EB8")
  ) +
  labs(
    x = "",
    y = "LISI"
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("Raw", "RPCA")),
    label = "p.signif"
  ) +
  ylim(0, max(lisi_long$LISI) * 1.5) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )
ggsave(
  file.path(fig_dir, "lisi_boxplot.pdf"),
  p,
  width = 1.5,
  height = 2
)
