

plot_vaf_density_by_group <- function(
    dat,
    group_var = "institution",
    x_axis_lab_spacing = 0.1,
    plot_title = NULL
) {
  gg <- ggplot(
    dat,
    aes(
      x = mut_vaf, 
      color = .data[[group_var]], 
      fill = .data[[group_var]]
    )
  ) + 
    geom_density(
      alpha = 0.2, 
      bounds = c(0,1), 
      trim = F,
      outline.type = "full"
    ) + 
    geom_vline(
      xintercept = 0.5, 
      color = "black", 
      linetype = "12"
    ) + 
    scale_x_continuous(
      breaks = seq(0, 1, by = x_axis_lab_spacing),
      labels = paste0(seq(0, 1, by = x_axis_lab_spacing) * 100, "%"),
      expand = c(0,0),
      limits = c(0,1)
    ) + 
    scale_y_continuous(expand = expansion(add = c(0,0), mult = c(0, 0.05))) + 
    labs(
      x = "Variant Allele Frequency",
      y = "Mutations (density)",
      title = plot_title
    ) + 
    khroma::scale_color_okabeito() + 
    khroma::scale_fill_okabeito() + 
    theme_classic() + 
    theme(
      plot.title.position = "plot",
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  return(gg)
  
}
