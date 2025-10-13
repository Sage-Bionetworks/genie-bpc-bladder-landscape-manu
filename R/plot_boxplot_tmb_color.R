plot_boxplot_tmb_color <- function(
  dat,
  x_var,
  x_lab = NULL,
  fill_var = "institution",
  plot_title = NULL,
  plot_subtitle = NULL,
  add_to_x = 0.1,
  hide_fill_legend = T
) {
  gg <- ggplot(
    data = dat,
    aes(
      x = .data[[x_var]] + add_to_x,
      y = panel_lab,
      fill = .data[[fill_var]]
    )
  ) +
    stat_boxplot(outlier.shape = NA, coef = 100) +
    # geom_jitter(width = 0, height = 1) +
    theme_classic() +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = x_lab
    ) +
    theme(
      plot.title.position = "plot",
      axis.title.y = element_blank(),
      axis.text.y = element_text(hjust = 0),
      strip.text = element_text(hjust = 0),
      legend.position = "bottom",
      axis.ticks.length.x = unit(0.2, 'cm')
    ) +
    scale_x_continuous(
      expand = expansion(mult = 0.01, add = 0),
      trans = 'log10',
      minor_breaks = log_tick_helper()
    ) +
    scale_y_discrete(limits = rev) +
    facet_wrap(
      vars(sample_type_simple_f),
      ncol = 1,
      scales = "free" # easiest way to get scales added in.
    ) +
    coord_cartesian(
      xlim = c(NA, max(dat[[x_var]] + add_to_x * 2))
    ) +
    scale_fill_vibrant() +
    guides(x = guide_axis(minor.ticks = T))

  if (hide_fill_legend) {
    gg <- gg + guides(fill = "none")
  }

  return(gg)
}
