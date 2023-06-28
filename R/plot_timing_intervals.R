# Old violin fill/line was cceeff and 225555
# Fill colors below are washed out versions of paul tol high contrast (the lines)

plot_timing_intervals <- function(
  dat, 
  time_var = "years",
  time_var_lab = "Years",
  facet_var = "interval",
  color_var = "origin",
  label_digits = 1,
  label_size = 2.5,
  label_stat_color = "black",
  violin_fill = c("#a9c6e2", "#f6e3aa", "#f5c6d0", "gray80"),
  violin_line = c("#004488", "#997700", "#994455", "gray20"),
  facet_scale = "fixed"
) {
  
  ld <- label_digits # more concise
  
  dft_sum <- dat %>%
    group_by(.data[[facet_var]]) %>%
    summarize(
      m = mean(.data[[time_var]], na.rm = T),
      std_dev = sd(.data[[time_var]], na.rm = T),
      med = median(.data[[time_var]], na.rm = T),
      fqt = quantile(.data[[time_var]], na.rm = T, probs = c(0.25)),
      tqt = quantile(.data[[time_var]], na.rm = T, probs = c(0.75)),
      miss_n = sum(is.na(.data[[time_var]]), na.rm = T),
      miss_prop = mean(is.na(.data[[time_var]]), na.rm = T),
      .groups = "drop"
    ) %>%
    mutate(
      str_mean = glue("Mean  (SD): {form_f(m, ld)} ({form_f(std_dev, ld)})"),
      str_med = glue("Median  (IQR): {form_f(med, ld)} ({form_f(fqt, ld)}, {form_f(tqt, ld)})"),
      str_miss = glue("Missing  (%): {miss_n} ({form_f(miss_prop*100, label_digits)}%)"),
      str_all = paste(str_mean, "<br>", str_med, "<br>", str_miss)
    )
  dft_sum <- select(dft_sum, all_of(facet_var), str_all)
  
  dat <- left_join(dat, dft_sum, by = facet_var) %>%
    # Remove all except the first row to avoid overplotting.
    group_by(.data[[facet_var]]) %>%
    mutate(str_all = case_when(
      row_number() %in% 1 ~ str_all,
      T ~ NA_character_
    ))
  
  x_max <- dat %>% pull(all_of(time_var)) %>% max(na.rm = T)
  
  gg <- ggplot(
    data = dat,
    aes(x = .data[[time_var]], y = 1)) + 
    geom_vline(xintercept = 0, color = "black", size = 0.5) + 
    geom_violin(
      aes(fill = .data[[color_var]], color = .data[[color_var]]),
      draw_quantiles = c(0.25, 0.5, 0.75)
    ) + 
    geom_jitter(width = 0, height = 0.25, alpha = 0.2, size = 0.25) + 
    geom_richtext(
      aes(x = x_max*.95, y = 1.05, label = str_all),
      vjust = 0, hjust = 1, size = label_size,
      color = label_stat_color,
      label.color = NA, label.padding = grid::unit(rep(0,4),"pt")
    ) + 
    facet_wrap(vars(interval), ncol = 1, scales = facet_scale) + 
    theme_bw() + 
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      strip.text.x = element_text(hjust = 0),
      legend.position = "none"
    ) + 
    scale_x_continuous(
      name = time_var_lab,
      n.breaks = 8, 
      expand = expansion(add = 0, mult = c(0, 0.05))
    ) +
    scale_color_manual(values = violin_line, drop = F) + 
    scale_fill_manual(values = violin_fill, drop = F)
  
  return(gg)
  
}
