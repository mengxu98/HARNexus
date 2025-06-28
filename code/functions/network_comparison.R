
reorder_data <- function(values) {
  names(values) <- methods_order
  values[methods_order]
}

create_split_data <- function(methods_order, values, break_points) {
  df <- data.frame(
    Method = factor(names(values), levels = methods_order),
    Value = unname(values)
  )

  df$type <- NA
  df$type[df$Value < break_points[1]] <- "small"
  df$type[df$Value > break_points[2]] <- "big"

  df_small <- filter(df, type == "big") %>%
    mutate(type = "small", Value = break_points[1]) %>%
    bind_rows(df)

  return(df_small)
}

create_split_plot <- function(
    data,
    title,
    y_lab,
    fill_color,
    break_points = NULL) {
  mymin <- function(y) ifelse(y <= break_points[1], 0, break_points[2])

  ggplot(data, aes(x = Method, y = Value, fill = Method)) +
    geom_rect(
      aes(
        xmin = as.numeric(Method) - 0.3,
        xmax = as.numeric(Method) + 0.3,
        ymin = mymin(Value),
        ymax = Value
      )
    ) +
    scale_fill_manual(values = fill_color) +
    geom_text(aes(label = format(Value, big.mark = ",")),
      vjust = 1, size = 3
    ) +
    facet_grid(type ~ ., scales = "free_y") +
    theme_bw() +
    theme(
      strip.text = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_y_continuous(
      labels = function(x) {
        format(x, scientific = FALSE, big.mark = ",", digits = 1)
      }
    ) +
    labs(title = title, y = y_lab, x = "")
}
