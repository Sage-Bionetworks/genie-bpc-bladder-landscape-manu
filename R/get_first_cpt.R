get_first_cpt <- function(
  ca_ind_dat,
  cpt_dat,
  include_sample_id = F,
  unit = "year"
) {
  rtn <- ca_ind_dat %>%
    select(record_id, ca_seq) %>%
    left_join(., cpt_dat, by = c("record_id", "ca_seq")) %>%
    arrange(cpt_number) %>%
    group_by(record_id) %>%
    slice(1) %>%
    ungroup()

  if (include_sample_id) {
    rtn_var <- c('record_id', 'ca_seq', 'cpt_genie_sample_id')
  } else {
    rtn_var <- c('record_id', 'ca_seq')
  }

  if (unit %in% 'year') {
    rtn %<>% select(all_of(rtn_var), dx_cpt_rep_yrs)
  } else if (unit %in% 'day') {
    rtn %<>% select(all_of(rtn_var), dx_cpt_rep_days)
  } else {
    cli_abort("Invalid unit, must be 'day', or 'year'")
  }

  return(rtn)
}
