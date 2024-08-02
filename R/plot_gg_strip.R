# A ggplot2 version of mice::stripplot(), which I wrote because the built in
#   version wasn't working well with categorical variables.
plot_gg_strip <- function(
    mids_obj,
    var,
    jitter = 0.25,
    pt_size = 2
) {
  dat <- mids_to_long_with_imputed(mids_obj = mids_obj, var = var)
  
  gg <- ggplot(
    data = dat,
    aes(x = .imp, y = .data[[var]], color = .imputed)
  ) + 
    geom_jitter(width = jitter, height = jitter, size = pt_size, alpha = 0.8) +
    scale_color_manual(values = c("#5285c9", "#C04C79")) + 
    theme_classic() +
    scale_x_continuous(breaks = 0:1000) + 
    labs(x = "Imputation iteration")
  
  return(gg)
}