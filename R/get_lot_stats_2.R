get_lot_stats_2 <- function(
    dat_reg,
    max_line = 3,
    time_vars = c("dx_reg_start_int_yrs", "dx_reg_end_all_int_yrs")
) {
  
  dat_reg %<>%
    filter(line_therapy <= max_line)
  
  dat_time <- dat_reg %>%
    group_by(line_therapy) %>%
    summarize(
      across(
        .cols = all_of(time_vars),
        .fns = list(
          q25 = (function(x) quantile(x, probs = 0.25, na.rm = T)),
          q50 = (function(x) quantile(x, probs = 0.50, na.rm = T)),
          q75 = (function(x) quantile(x, probs = 0.75, na.rm = T))
        )
      ),
      .groups = "drop"
    )
  
  dat_overall <- dat_reg %>%
    group_by(line_therapy) %>%
    summarize(n = n())
  
  dat_rtn <- full_join(
    dat_overall,
    dat_time,
    by = "line_therapy"
  )
  
  dat_rtn %<>%
    mutate(line_therapy_txt = english::ordinal(line_therapy)) %>%
    select(-line_therapy) %>%
    select(line_therapy_txt, everything())
  
  return(dat_rtn)
}
