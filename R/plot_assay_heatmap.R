plot_assay_heatmap <- function(
    dat,
    plot_title = "Gene panel coverage" 
) {
  gg <- ggplot(
    data = dat,
    aes(x = hugo, y = assay_lab)
  ) + 
    geom_raster() +
    theme_classic() + 
    scale_y_discrete(limits = rev) + 
    labs(x = "Gene (one line per gene)",
         title = plot_title,
         subtitle = "Panels annotated with [number of genes] and (patients, samples) tested.") + 
    theme(
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      plot.title.position = "plot"
    )
  
  return(gg)
  
}