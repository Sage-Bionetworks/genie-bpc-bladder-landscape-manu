forest_mod_natural_scale <- function(
    dat,
    model_var = "model"
) {
  
  gg <- ggplot(
    dat,
    aes(x = estimate, xmin = conf.low, xmax = conf.high, y = term,
        color = .data[[model_var]])
  ) + 
    geom_vline(color = 'gray70', linewidth = 2, alpha = 0.5, xintercept = 0) + 
    geom_pointrange(position = position_dodge2(width = 0.5),
                    shape = 124) + 
    theme_bw() + 
    scale_color_highcontrast() +
    guides(color = guide_legend(title = NULL)) +
    theme(legend.position = "bottom")
  
  return(gg)
  
}