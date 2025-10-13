# Adds the cohort entry time, which is when both the index event and the
#   truncation time variable have both happened.
# In GENIE a typical index event would be regimen start time.
# Truncation time variable is the first cancer panel test.
# This function assumes the index interval is in days, and from diagnosis (unchecked)
add_entry_time <- function(
  dat_index,
  var_index,
  dat_ca_ind,
  dat_cpt
) {
  dft_first_cpt <- get_first_cpt(
    dat_ca_ind,
    dat_cpt,
    unit = 'day',
    include_sample_id = F
  ) %>%
    rename(dx_first_cpt_rep_days = dx_cpt_rep_days)

  rtn <- left_join(
    dat_index,
    dft_first_cpt,
    by = c("record_id", 'ca_seq')
  )

  rtn %<>%
    mutate(
      # risk set entry happens when BOTH they have a CPT test and take a post-met drug.
      dx_entry = case_when(
        # they have to HAVE both those events, otherwise they did not enter:
        is.na(dx_first_cpt_rep_days) | is.na(.data[[var_index]]) ~ NA_real_,
        T ~ pmax(dx_first_cpt_rep_days, .data[[var_index]], na.rm = T)
      )
    )

  return(rtn)
}
