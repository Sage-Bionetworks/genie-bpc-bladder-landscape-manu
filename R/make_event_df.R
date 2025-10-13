# cohort_dat describes the people who should be included.
make_event_df <- function(ca_ind_dat, cpt_dat, cohort_dat = ca_ind_dat) {
  # for each index cancer get the first cpt
  first_panel_dates <- get_first_cpt(ca_ind_dat, cpt_dat) %>%
    select(record_id, t_yrs = dx_cpt_rep_yrs)

  event_dat <- cohort_dat %>%
    select(record_id, ca_seq) %>%
    left_join(
      .,
      select(ca_ind_dat, record_id, ca_seq, tt_os_dx_yrs, os_dx_status),
      by = c("record_id", "ca_seq")
    )

  event_dat %<>%
    mutate(
      event = case_when(
        os_dx_status %in% 0 ~ "Censored",
        os_dx_status %in% 1 ~ "Death"
      )
    ) %>%
    select(-os_dx_status) %>%
    rename(t_yrs = tt_os_dx_yrs)

  event_dat %<>%
    bind_rows(., first_panel_dates) %>%
    mutate(
      event = if_else(is.na(event), "First NGS", event),
      event = factor(event, levels = c("First NGS", "Censored", "Death"))
    )

  return(event_dat)
}
