get_quantile_df <- function(dat, var, probs = c(0.25, 0.5, 0.75)) {
  y <- quantile(dat[[var]], probs = probs)
  rtn <- tibble(y = y) %>%
    mutate({{ var }} := probs)

  return(rtn)
}
