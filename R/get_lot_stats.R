get_lot_stats <- function(
    dat_reg,
    max_reg_number = 3,
    time_vars = c("dx_reg_start_int_yrs", "dx_reg_end_all_int_yrs")
) {
  
  dat_reg <- lot_filter_prep(dat_reg, max_reg = max_reg_number)
  
  dat_time <- dat_reg %>%
    group_by(.conseq_reg_number) %>%
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
    group_by(.conseq_reg_number) %>%
    summarize(n = n())
  
  dat_rtn <- full_join(
    dat_overall,
    dat_time,
    by = ".conseq_reg_number"
  )
  
  dat_rtn %<>%
    mutate(line_of_therapy = english::ordinal(.conseq_reg_number)) %>%
    select(-.conseq_reg_number) %>%
    select(line_of_therapy, everything())
  
  return(dat_rtn)
}
