# ca_ind_dat is the index cancer diagnosis dataset.
make_dmet_musc_prop_status_block <- function(
    ca_ind_dat,
    truncate_musc_inv_at_dx = T
) {
  
  ci_sub <- ca_ind_dat %>% 
    # just to be clear about variables we're using:
    select(
      record_id, 
      ca_seq,
      stage_dx_iv,
      ca_dmets_yn,
      dx_to_dmets_yrs,
      dx_muscprop_invasive_yrs,
      tt_os_dx_yrs
    )
  
  last_fu_time <- ci_sub %>%
    # loss to follow up:
    mutate(event = "ltfu") %>%
    rename(dx_to_event_yrs = tt_os_dx_yrs) %>%
    select(record_id, ca_seq, event, dx_to_event_yrs)
  
  dmet_at_onset <- ci_sub %>%
    filter(stage_dx_iv %in% "Stage IV" & ca_dmets_yn %in% "Yes") %>%
    mutate(event = "dmet") %>%
    mutate(dx_to_event_yrs = 0) %>%
    select(record_id, ca_seq, event, dx_to_event_yrs)
  
  dmet_in_fu <- ci_sub %>%
    filter(!(stage_dx_iv %in% "Stage IV" & ca_dmets_yn %in% "Yes")) %>%
    filter(!is.na(dx_to_dmets_yrs)) %>%
    mutate(event = "dmet") %>%
    rename(dx_to_event_yrs = dx_to_dmets_yrs) %>%
    select(record_id, ca_seq, event, dx_to_event_yrs)
  
  musc_inv <- ci_sub %>%
    select(record_id, ca_seq, dx_muscprop_invasive_yrs) %>%
    filter(!is.na(dx_muscprop_invasive_yrs)) %>%
    mutate(event = "musc_inv") %>%
    rename(dx_to_event_yrs = dx_muscprop_invasive_yrs) %>%
    select(record_id, ca_seq, event, dx_to_event_yrs) 
  
  if (truncate_musc_inv_at_dx) {
    musc_inv %<>%
      mutate(
        dx_to_event_yrs = if_else(dx_to_event_yrs < 0, 0, dx_to_event_yrs)
      )
  }
  
  dx <- ci_sub %>%
    select(record_id, ca_seq) %>%
    mutate(event = "dx", dx_to_event_yrs = 0)
  
  event_priority <- c("dx", "musc_inv", "dmet", "ltfu")
  
  event_df <- bind_rows(
    last_fu_time,
    dmet_at_onset,
    dmet_in_fu,
    musc_inv,
    dx
  ) %>%
    mutate(event = factor(event, levels = event_priority)) %>%
    arrange(record_id, ca_seq, dx_to_event_yrs, event) 
  
  # if there are multiple events at the exact same time, take the highest
  #   priority event.
  event_df %<>%
    group_by(record_id, ca_seq, dx_to_event_yrs) %>%
    slice(n()) %>%
    ungroup(.)
  
  event_df %>%
    filter(record_id %in% 'GENIE-DFCI-000095' & ca_seq %in% 3) %>%
    print(.)
  
  event_df %<>%
    mutate(
      dx_block_start = dx_to_event_yrs,
      dx_block_end = dplyr::lead(dx_to_event_yrs)
    )
  
  event_df %>%
    filter(record_id %in% 'GENIE-DFCI-000095' & ca_seq %in% 3) %>%
    print(.)
  
  # The last block should be just the end of followup time.  We've already
  #   used that information with the lead() call above.
  event_df %<>%
    group_by(record_id, ca_seq) %>%
    slice(1:(n()-1)) %>%
    ungroup(.)
  
  event_to_block <- tribble(
    ~event, ~status,
    "dx", "Non-invasive",
    "musc_inv", "Invasive", 
    "dmet", "Metastatic"
  )
  
  rtn <- event_df %>%
    left_join(., event_to_block, by = "event") %>%
    mutate(status = factor(status, levels = event_to_block$status)) %>%
    select(record_id, ca_seq, status, dx_block_start, dx_block_end)
  
  rtn %>%
    filter(record_id %in% 'GENIE-DFCI-000095' & ca_seq %in% 3) %>%
    print(.)
      
  return(rtn)
  
  
  
}


