source("code/functions/prepare_env.R")

fig_dir <- check_dir("figures/har_tf")

tfs <- read.csv("results/har_tf/tfs.csv")[, 1]
res <- RunEnrichment(
  geneID = tfs,
  species = "Homo_sapiens",
  db = "GO_BP",
  Ensembl_version = 113
)


p1 <- EnrichmentPlot(
  res = res,
  plot_type = "lollipop"
)
ggsave(
  file.path(fig_dir, "tfs_enrichment_lollipop.pdf"),
  p1,
  width = 7.5,
  height = 4
)

EnrichmentPlot(
  res = res,
  plot_type = "network"
)
ggsave(
  file.path(fig_dir, "tfs_enrichment_network.pdf"),
  width = 12,
  height = 10
)

EnrichmentPlot(
  res = res,
  plot_type = "enrichmap",
  theme_use = "theme_blank",
  enrlichmap_nlabel = 5,
  theme_args = list(
    plot.subtitle = ggplot2::element_text(size = 10),
    strip.text = ggplot2::element_text(size = 8)
  )
)
ggsave(
  file.path(fig_dir, "tfs_enrichment_enrichmap_blank.pdf"),
  width = 12,
  height = 10
)
