plot_ecdf_by_panel <- function(
    dat, 
    x_var,
    x_lab = NULL,
    plot_title = NULL,
    plot_subtitle = NULL,
    add_to_x = 0.1,
    x_int_lines = c(10.1)
) { 
  
  gg <- ggplot(
    data = dat,
    aes(
      x = .data[[x_var]] + add_to_x, 
      color = institution
    )
  ) + 
    stat_ecdf() + 
    geom_vline(xintercept = x_int_lines) + 
    theme_bw() + 
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = x_lab
    ) +
    theme(
      plot.title.position = "plot",
      axis.title.y = element_blank(),
      axis.text.y = element_text(hjust = 0),
      strip.text = element_text(hjust = 0)
    ) + 
    scale_x_continuous(
      expand = expansion(mult = 0.01, add = 0),
      trans = 'log10'
    ) +
    scale_y_continuous(
      breaks = seq(0,1, by = 0.5),
      labels = paste0(seq(0,1, by = 0.5)*100, "%"),
      minor_breaks = seq(0,1, by = 0.25)
    ) + 
    facet_wrap(
      vars(sample_type_simple_f), 
      ncol = 1,
      scales = "free" # free_x is the easiest way to get axes.
    ) +
    coord_cartesian(
      xlim = c(NA, max(dat[[x_var]] + add_to_x * 2))
    ) + 
    scale_color_vibrant()
  
  return(gg)
  
}
