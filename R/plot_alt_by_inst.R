
plot_alt_by_inst <- function(
    dat,
    prop_var = "prop_alt",
    plot_title = "Alteration frequencies by institution",
    plot_subtitle = "Any alteration - oncogenic or not.",
    reorder_by_prop = T,
    point_size = 1
) {
  
  dat %<>% mutate(
    shape_scale_help = factor(
      institution %in% "All",
      levels = c(F, T)
    )
  )
  
  if (reorder_by_prop) {
    # only works because "all" is first:
    dat %<>%
      arrange(institution, desc(.data[[prop_var]])) %>%
      mutate(hugo = fct_inorder(hugo))
  }
  
  gg <- ggplot(
    data = dat,
    aes(
      x = .data[[prop_var]],
      y = hugo,
      color = institution,
      shape = shape_scale_help,
      size = shape_scale_help
    )
  ) + 
    geom_hline(yintercept = 0:100, color = "gray90") + 
    geom_point(alpha = 0.8) + 
    scale_x_continuous(
      breaks = seq(0, 100, by = 10)/100,
      minor_breaks = seq(0, 100, by = 5)/100,
      labels = paste(seq(0, 100, by = 10), "%"),
      name = "Proportion altered",
      expand = expansion(add = 0, mult = 0.02)
    ) + 
    scale_shape_manual(
      values = c(19, 15),
      guide = 'none'
    ) + 
    scale_color_manual(
      values = c("gray20", '#ee7733', '#0077bb', '#33bbee', '#ee3377')
    ) + 
    # point_size is a scaling factor here.
    scale_size_manual(values = point_size*c(1,2), guide = "none") + 
    scale_y_discrete(
      limits = rev,
      name = NULL
    ) + 
    labs(
      title = plot_title,
      subtitle = plot_subtitle
    ) + 
    theme_bw() + 
    theme(
      legend.position = 'bottom',
      plot.title.position = "plot"
    )
  
  return(gg)
  
}
