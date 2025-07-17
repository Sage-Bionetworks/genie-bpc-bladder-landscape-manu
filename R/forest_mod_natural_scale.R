forest_mod_natural_scale <- function(
  dat,
  model_var = "model",
  term_var = "term",
  pal = NULL
) {
  gg <- ggplot(
    dat,
    aes(
      x = estimate,
      xmin = conf.low,
      xmax = conf.high,
      y = .data[[term_var]],
      color = .data[[model_var]]
    )
  ) +
    geom_vline(color = 'gray20', linewidth = 0.5, alpha = 1, xintercept = 0) +
    geom_pointrange(position = position_dodge2(width = 0.5), shape = 124) +
    theme_bw() +
    guides(color = guide_legend(title = NULL)) +
    theme(legend.position = "bottom", axis.text.y = element_markdown())

  if (is.null(pal)) {
    gg <- gg +
      scale_color_bmj()
  } else {
    gg <- gg +
      scale_color_manual(values = pal)
  }

  return(gg)
}
