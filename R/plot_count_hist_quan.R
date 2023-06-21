#' @title Plot observation count histograms with quantiles
#' 
#' @param dat A dataframe with required columns type_f and n.
#' @param hist_col Color of the histogram bars.
#' @param quantiles_to_plot Quantiles to put vertical bars under.
#' @param pal_quantiles Colors of the quantile vertical bars.
#' @param x_lab Label of the X axis (default is "Number of Measurements").
#' @param y_lab Label of the Y axis (default is "Subjects").
#' @param title Label for the title.  Defaults to quantile label for default quantiles.  If no title is desired, use `title = NULL`
#' @param y_bottom_exp Expansion of the bottom of the Y axis (multiplicative).
plot_count_hist_quan <- function(
  dat, 
  hist_col = "gray20", 
  quantiles_to_plot = c(0.25, 0.5, 0.75),
  quantile_names = c("1st quartile", "Median", "3rd quartile"),
  pal_quantiles = c('#0077bb', '#ee7733', '#ee3377'),
  plot_title = "default",
  x_lab = "Number of Measurements",
  y_lab = "Subjects",
  y_bottom_exp = 0.05,
  plot_title_size = 11
) {
  
  if (length(quantiles_to_plot) > length(pal_quantiles)) {
    cli::cli_abort("Length of quantiles_to_plot exceeds pal_quantiles.")
  } else if (length(pal_quantiles) > length(quantiles_to_plot)) {
    pal_quantiles <- pal_quantiles[1:length(quantiles_to_plot)]
  }
  
  if (plot_title %in% "default") {
    plot_title <- tibble(
      name = quantile_names, 
      val = pal_quantiles
    ) %>% 
      mutate(
        name = glue("<span style='color:{val}'>{name}</span>")
      ) %>%
      pull(name) %>%
      paste(collapse = ", ") %>%
      paste0("Histograms annotated with ", ., ".")
  }
  
  
  
  gg <- ggplot(data = dat,
               aes(x = n)) + 
    # this part does the quantile summaries
    stat_summary(
      geom = "vline",
      orientation = "y",
      # y is a required aesthetic, so use a dummy value
      aes(y = 1, xintercept = after_stat(x)),
      fun = function(x) {
        quantile(x, probs = c(0.25, 0.5, 0.75))
      },
      alpha = 1, linewidth = 0.5,
      color = rep(
        pal_quantiles,
        times = length(unique(dat$type_f))
      )
    ) + 
    
    geom_histogram(binwidth = 1, color = hist_col, fill = hist_col) + 
    theme_bw() +
    labs(title = plot_title) + 
    facet_wrap(vars(type_f), scales = "free", nrow = 2) + 
    theme(
      strip.text = element_text(hjust =0),
      panel.grid.minor.x = element_blank(),
      axis.title.y = element_markdown(angle = 0, vjust = 0.5),
      plot.title = element_markdown(size = plot_title_size),
      plot.title.position = "plot"
    ) + 
    scale_x_continuous(
      name = x_lab,
      expand = expansion(add = 0, mult = 0),
      n.breaks = 6
    ) +
    scale_y_continuous(
      name = y_lab,
      expand = expansion(add = c(0,0), mult = c(y_bottom_exp,0.05))
    ) 
  return(gg)
  
}