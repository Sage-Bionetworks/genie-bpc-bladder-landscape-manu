get_binary_feature_pos <- function(
    dat,
    ignore_cols = character(0)
) {
  dat %<>%
    select(
      -all_of(ignore_cols)
    ) 
  
  chk_binary <- dat %>%
    summarize(
      across(
        .cols = everything(),
        .fns = ~mean(.x %in% c(0,1) | is.na(.x))
      )
    ) %>%
    pivot_longer(
      cols = everything(),
      names_to = "feature",
      values_to = "prop_binary_or_na"
    )
  
  if (any(chk_binary$prop_binary_or_na < 1)) {
    str_bin_error <- chk_binary %>% 
      filter(prop_binary_or_na < 1) %>%
      arrange(prop_binary_or_na) %>%
      slice(min(3,n())) %>%
      pull(feature) %>%
      paste(collapse = ", ")
    cli_abort("Not all features are binary or NA: {str_bin_error}")
  }
  
  rtn <- dat %>%
    summarize(
      across(
        .cols = everything(),
        .fns = ~sum(.x %in% 1)
      )
    ) %>%
    pivot_longer(
      cols = everything(),
      names_to = "feature",
      values_to = "num_pos"
    ) %>%
    arrange(desc(num_pos))
  
  return(rtn)
}