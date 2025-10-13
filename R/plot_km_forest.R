plot_km_forest <- function(
  dat,
  y,
  plot_infinite = T
) {
  if (plot_infinite) {
    # take both the estimate and high confidence interval in case there is an
    #   estimate that exceeds all the high confidence intervals.
    ci_max <- max(c(dat$conf.high, dat$median), na.rm = T)
    inf_dat <- dat %>%
      filter(is.na(conf.high) & !is.na(conf.low)) %>%
      mutate(conf.high = ci_max * 1.05)

    dat <- dat %>%
      mutate(
        est_type = if_else(
          is.na(conf.high) | is.na(conf.low),
          "infinite",
          "finite"
        )
      )

    inf_dat <- dat %>%
      filter(est_type %in% "infinite" & !is.na(conf.low)) %>%
      mutate(conf.high = ci_max * 1.05)

    dat %<>%
      replace_na(list(conf.high = ci_max * 1.05))
  }

  gg <- ggplot(
    dat,
    aes(
      xmin = conf.low,
      xmax = conf.high,
      x = median,
      y = .data[[y]],
      color = est_type
    )
  ) +
    theme_bw() +
    geom_pointrange() +
    geom_point(data = inf_dat, shape = 1, aes(x = conf.high, y = .data[[y]])) +
    scale_color_jama() +
    theme(
      plot.title.position = 'plot',
      legend.position = 'bottom'
    )

  return(gg)
}
