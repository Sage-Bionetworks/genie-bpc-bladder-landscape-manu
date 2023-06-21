# ca_ind_dat is the index cancer diagnosis dataset.
make_dmet_status_block <- function(ca_ind_dat) {
  
  dmet_stat_levs <- c("No dmet noted",
                      "Distant Metastasis")
  
  ci_sub <- ca_ind_dat %>% 
    # just to be clear about variables we're using:
    select(
      record_id, 
      ca_seq,
      stage_dx_iv,
      ca_dmets_yn,
      dx_to_dmets_yrs,
      tt_os_dx_yrs
    )
  
  dmet_at_onset <- ci_sub %>%
    filter(stage_dx_iv %in% "Stage IV" & ca_dmets_yn %in% "Yes") %>%
    mutate(dmet_status = dmet_stat_levs[2],
           dx_block_start = 0,
           dx_block_end = tt_os_dx_yrs) 
  
  dmet_in_fu <- ci_sub %>%
    # A note on this filtering: There is at least one person who was diagnosed
    #   as stage 4 but has dmets_stage_i_iii == 1.  Instead of messing with
    #   those variables, we're just taking the complement of the previous filter
    #   and then filtering down to anyone with a dmet time listed - that way
    #   we get everyone with usable information.
    filter(!(stage_dx_iv %in% "Stage IV" & ca_dmets_yn %in% "Yes")) %>%
    filter(!is.na(dx_to_dmets_yrs)) %>%
    # for each person in this subset we need 2 blocks: dx to dmet and dmet to death/censor.
    group_by(record_id, ca_seq) %>%
    slice(rep(1:n(), each = 2)) %>%
    mutate(
      dmet_status = if_else(row_number() %in% 1,
                            dmet_stat_levs[1],
                            dmet_stat_levs[2]),
      dx_block_start = if_else(row_number() %in% 1, 0, dx_to_dmets_yrs),
      dx_block_end = if_else(row_number() %in% 1, dx_to_dmets_yrs, tt_os_dx_yrs)
    ) 
  
  dmet_never <- ci_sub %>%
    # See note on filtering above.
    filter(!(stage_dx_iv %in% "Stage IV" & ca_dmets_yn %in% "Yes")) %>%
    filter(is.na(dx_to_dmets_yrs)) %>%
    mutate(dmet_status = dmet_stat_levs[1],
           dx_block_start = 0,
           dx_block_end = tt_os_dx_yrs) 
  
  dmet_all <- bind_rows(
    dmet_never,
    dmet_in_fu,
    dmet_at_onset
  ) %>%
    # keep tt_os_dx_yrs for sorting the plot.
    select(record_id, ca_seq, dmet_status, contains("dx_block"), tt_os_dx_yrs) %>%
    mutate(dmet_status = factor(dmet_status, dmet_stat_levs))
  
  return(dmet_all)
  
}