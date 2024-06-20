#
#' @param drug_dat A dataframe with columns "drug" and "n".  Used in the order given.
plot_drug_prop <- function(drug_dat, cohort_n, pal = NULL, pal_seed = 278, plot_title = NULL, x_expansion = 0) {
  if (is.null(pal)) {
    set.seed(pal_seed)
    pal <- sample(make_sun_pal(nrow(drug_dat)))
  }
  
  drug_dat <- drug_dat %>%
    select(drug, n) %>%
    mutate(
      prop = n/cohort_n,
      str = glue("n={n} ({formatC(prop*100, digits = 1, format = 'f')}%)")
    ) %>%
    mutate(drug = forcats::fct_inorder(drug)) 
  
  gg <- ggplot(data = drug_dat,
               aes(x = prop, y = drug, fill = drug)) +
    geom_col() + 
    annotate(
      geom = 'rect', 
      xmin = 1, xmax = 1.5, ymin = -2, ymax = 15,
      color = "gray80", fill = "gray80"
    ) + 
    geom_text(size = 2.5, 
              nudge_x = .05, 
              hjust = 0, 
              aes(label = str)) + 
    labs(title = plot_title) + 
    scale_fill_manual(
      values = pal
    ) + 
    scale_x_continuous(
      name = "Proportion of cohort",
      breaks = c(0,0.5,1),
      labels = paste(c(0, 50, 100), "%"),
      expand = expansion(mult = c(0,0), add = c(0,0))
    ) + 
    scale_y_discrete(limits = rev) + 
    coord_cartesian(xlim = c(0,1 + x_expansion), ylim = c(0.9, 10.1)) + 
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      legend.position = "none",
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.title.position = 'plot',
      # Allows the "%" to stay on "100%"
      plot.margin = margin(r = 20, unit = "pt")
    )
  
  return(gg)
}





#' @param drug_dat A dataframe with columns "drug" and "n".  Used in the order given.
plot_drug_prop_vert <- function(
    drug_dat, 
    cohort_n, 
    pal = NULL,
    pal_seed = 278, 
    plot_title = NULL,
    x_expansion = 0,
    y_left_margin_exp = 20
) {
  if (is.null(pal)) {
    set.seed(pal_seed)
    pal <- sample(make_sun_pal(nrow(drug_dat)))
  }
  
  drug_dat <- drug_dat %>%
    select(drug, n) %>%
    mutate(
      prop = n/cohort_n,
    ) %>%
    mutate(drug = forcats::fct_inorder(drug)) 
  
  gg <- ggplot(data = drug_dat,
               aes(y = prop, x = drug, fill = drug)) +
    geom_col() + 
    labs(title = plot_title) + 
    scale_fill_manual(
      values = pal
    ) + 
    scale_x_discrete(
      name = NULL
    ) + 
    scale_y_continuous(
      name = "Proportion of cohort",
      n.breaks = 6,
      labels = scales::label_percent(),
      expand = expansion(mult = c(0,0.05), add = c(0,0))
    ) + 
#     coord_cartesian(xlim = c(0,1 + x_expansion), ylim = c(0.9, 10.1)) + 
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.title.position = 'plot',
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.margin = margin(l = y_left_margin_exp, unit = "pt") # for Gemcitabine.
    )
  
  return(gg)
}










