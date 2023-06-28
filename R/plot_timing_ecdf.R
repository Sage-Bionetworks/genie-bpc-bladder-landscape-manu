


plot_timing_ecdf <- function(
  dat, 
  time_var = "years",
  time_var_lab = "Years",
  facet_var = "interval",
  color_var = "origin",
  pal_line = c("#004488", "#997700", "#994455", "gray20"),
  facet_scale = "fixed"
) {
  
  gg <- ggplot(
    data = dat,
    aes(x = .data[[time_var]], y = 1)) + 
    geom_vline(xintercept = 0, color = "gray70", size = 0.5) + 
    stat_ecdf(
      pad = F,
      aes(color = .data[[color_var]])
    ) + 
    facet_wrap(vars(interval), ncol = 1, scales = facet_scale) + 
    theme_bw() + 
    theme(
      strip.text.x = element_text(hjust = 0),
      legend.position = "none"
    ) + 
    scale_x_continuous(
      name = time_var_lab,
      n.breaks = 8, 
      expand = expansion(add = 0, mult = c(0, 0.05))
    ) +
    scale_y_continuous(
      position = "right", name = NULL,
      breaks = seq(0,1,by = 0.25),
      labels = c("0%", "", "50%", "", "100%"),
      expand = expansion(add = 0, mult = 0.01)
    ) + 
    scale_color_manual(values = pal_line)
  
  return(gg)
}