
test_kendall_tau <- function(
    dat,
    var_trunc,
    var_event_time,
    var_event_ind, # could be modified to work without this if desired.
    method = "MB"
) {
  
  res <- tranSurv::cKendall(
    trun = dat[[var_trunc]],
    obs = dat[[var_event_time]],
    delta = dat[[var_event_ind]],
    method = method
  ) 
  
  res %<>% 
    `[`(., c("PE", "SE", "STAT", "p.value", "trans", "a")) %>%
    as_tibble(.) %>%
    mutate(
      var_trunc = var_trunc,
      var_event_time = var_event_time,
      method = method
    ) %>%
    select(
      var_trunc:method,
      # adopting broom's standard stuff.
      estimate = "PE",
      std.error = "SE",
      statistic = "STAT",
      p.value
    ) 
  
  
  return(res)
  
}
