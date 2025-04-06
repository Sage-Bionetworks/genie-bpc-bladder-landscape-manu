univariate_logistic_help <- function(dat, y, x) {
  form_str = paste0(y, "~", x)
  
  res <- glm(
    data = ddr_outcome,
    formula = as.formula(form_str),
    family = 'binomial'
  ) %>%
    broom::tidy(., conf.int = T) %>%
    filter(!(term %in% "(Intercept)"))
  
  return(res)
  
}