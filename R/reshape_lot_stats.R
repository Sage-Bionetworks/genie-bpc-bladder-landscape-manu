
reshape_lot_stats <- function(
    dat_lot,
    n_denom = NULL,
    digits_t = 2
) {
  rtn <- dat_lot %>%
    pivot_longer(
      cols = -line_therapy_txt
    ) %>%
    separate(
      col = name,
      into = c("name", "quantile"),
      sep = "_q"
    ) %>%
    mutate(
      quantile = paste0("q", quantile)
    ) %>%
    pivot_wider(
      names_from = quantile,
      values_from = value
    ) 
  
  rtn %<>% mutate(
    str = case_when(
      is.na(q25) ~ as.character(qNA),
      T ~ as.character(
        glue("{form_f(q50, digits_t)} ({form_f(q25, digits_t)}, {form_f(q75, digits_t)})")
      )
    )
  )
  
  rtn %<>%
    select(line_therapy_txt, name, str)
  
  if (!is.null(n_denom)) {
    # Show the percentage relative to the denominator if true.
    rtn %<>%
      mutate(
        str = if_else(
          name %in% "n",
          n_pct_str(as.numeric(str), n_denom, show_d = T, digits = 0),
          str
        )
      )
  }
  
  
  rtn %<>%
    pivot_wider(
      names_from = line_therapy_txt,
      values_from = str
    ) %>%
    rename(quantity = name) %>% 
    rename_all(
      stringr::str_to_title
    )
  
  return(rtn)
}
