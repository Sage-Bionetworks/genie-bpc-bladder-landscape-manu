beautify_strings <- function(var, dict, type = "default", factorize = T) {
  new_vec <- tibble(raw = var) %>%
    left_join(., dict, by = c("raw")) %>%
    pull(.data[[type]])
  
  if (factorize) {
    levs <- dict[[type]]
    new_vec <- factor(new_vec, levels = levs)
    new_vec <- forcats::fct_drop(new_vec)
  }
  
  return(new_vec)
}
