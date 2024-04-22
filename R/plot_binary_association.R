plot_binary_association <- function(
    dat,
    x_var,
    y_var,
    pval_var = "p.value",
    pval_thresh = 10^(-5),
    colorbar_width = 20,
    show_p_sig = F,
    alpha = 0.05,
    label_var = NULL,
    nudge_y = 0,
    nudge_x_lab = 0.15,
    label_size = 2
) {
  
  # for any p value under the threshold, raise it to the threshold.
  # The point of this is to avoid swamping the plot with crazy values like
  #  10^(-180) which will skew the color scale to insanity.
  dat %<>%
    mutate(
      {{pval_var}} := if_else(
        pval_thresh > .data[[pval_var]],
        pval_thresh,
        .data[[pval_var]]
      )
    )
  
  dat %<>%
    mutate(
      log10_pval = log10(.data[[pval_var]]),
      log10_pval_signed = if_else(
        odds_ratio > 1,
        -log10_pval, # conflicing (negative association, "mutually exclusive")
        log10_pval   # co-occuring (positive association)
      ),
      sig = if_else(.data[[pval_var]] < alpha, T, F, F)
    )
  
  biggest_log10_p <- max(abs(dat$log10_pval))
  
  gg <- ggplot(
    data = dat,
    aes(
      x = .data[[x_var]],
      y = .data[[y_var]],
      fill = log10_pval_signed,
      color = (if (show_p_sig) sig else NULL)
    )
  ) + 
    geom_tile(linewidth = 0) +
    geom_point(position = position_nudge(y = nudge_y)) + 
    theme_bw() + 
    theme(
      legend.position = "bottom",
      panel.grid = element_blank()
    ) + 
    scale_fill_continuous_diverging(
      h1 = 195, h2 = 325,
      c1 = 90,
      l1 = 50, l2 = 95,
      p1 = 1, p2 = 1,
      mid = 0,
      limits = c(-1,1) * ceiling(biggest_log10_p),
      breaks = c(-1,1) * ceiling(biggest_log10_p),
      labels = c(
        'Negatively Associated \n ("mututally exclusive")', 
        "Positively Associated \n (co-occuring)"
      ),
      name = NULL
    ) + 
    scale_x_discrete(
      drop = F, 
      position = "top",
      name = NULL
    ) + 
    scale_y_discrete(
      drop = F, 
      name = NULL
    ) +
    scale_color_manual(
      values = c("transparent", "black")
    ) +
    guides(
      fill = guide_colorbar(
        barwidth = colorbar_width
      ),
      color = "none"
    )
  
  if (!is.null(label_var)) {
    gg <- gg +
      geom_text(
        inherit.aes = F,
        aes(
          x = .data[[x_var]],
          y = .data[[y_var]],
          label = .data[[label_var]]
        ),
        color = "black",
        size = label_size,
        nudge_y = nudge_y,
        nudge_x = -0.5 + nudge_x_lab,
        hjust = 0,
        family = "mono"
      )
  }
  
  return(gg)
  
  
  
  
  
  
  
} 
