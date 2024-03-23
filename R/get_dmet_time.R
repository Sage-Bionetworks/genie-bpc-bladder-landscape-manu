# Get the dmet time for all subjects.
get_dmet_time <- function(
    ca_ind_dat,
    annotate_type = F
) {
  
  ci_sub <- ca_ind_dat %>% 
    # just to be clear about variables we're using:
    select(
      record_id, 
      ca_seq,
      stage_dx_iv,
      ca_dmets_yn,
      dmets_stage_i_iii,
      matches('^dx_to_dmets.*_yrs$')
    )
  
  dmet_at_onset <- ci_sub %>%
    filter(stage_dx_iv %in% "Stage IV" & ca_dmets_yn %in% "Yes") %>%
    mutate(dx_dmet_yrs = 0) %>%
    select(record_id, ca_seq, dx_dmet_yrs)
  
  dmet_later_stage_i_iii <- ci_sub %>%
    filter(!(stage_dx_iv %in% "Stage IV") & dmets_stage_i_iii %in% 1)
  
  chk_1 <- dmet_later_stage_i_iii %>%
    filter(xor(dmets_stage_i_iii %in% 1, !is.na(dx_to_dmets_yrs)))
  if (nrow(chk_1) > 0) {
    cli::cli_alert_danger("dmets_stage_i_iii and dx_to_dmets_yrs have a discrepancy")
  }
  
  dmet_later_stage_i_iii %<>%
    select(record_id, ca_seq, dx_dmet_yrs = dx_to_dmets_yrs)
  
  
  dmet_later_stage_iv <- ci_sub %>%
    filter(stage_dx_iv %in% "Stage IV" & !(ca_dmets_yn %in% "Yes")) %>%
    collapse_site_mets(., rtn_col_name = "dx_dmet_yrs") %>%
    filter(!is.na(dx_dmet_yrs))
  
  if (annotate_type) {
    dmet_at_onset %<>% mutate(.met_type = "stage_iv_with_immediate_dmet")
    dmet_later_stage_iv %<>% mutate(.met_type = "stage_iv_dmet_later")
    dmet_later_stage_i_iii %<>% mutate(.met_type = "not_stage_iv_dmet_later")
  }
  
  
  rtn <- bind_rows(
    dmet_at_onset,
    dmet_later_stage_i_iii,
    dmet_later_stage_iv
  )
  
  if (max(pull(count(rtn, record_id, ca_seq, sort = T), n), na.rm = T) > 1) {
    cli::cli_abort("Duplicated rows - clear error")
  }
  
  return(rtn)
  
  
  
}

# get_dmet_time(dft_ca_ind, annotate_type = T)



