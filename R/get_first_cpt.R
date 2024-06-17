get_first_cpt <- function(ca_ind_dat, cpt_dat, include_sample_id = F) {
  rtn <- ca_ind_dat %>% 
    select(record_id, ca_seq) %>%
    left_join(., cpt_dat, by = c("record_id", "ca_seq")) %>%
    arrange(cpt_number) %>%
    group_by(record_id) %>%
    slice(1) %>%
    ungroup() 
  
  if (include_sample_id) {
    rtn %<>% select(record_id, ca_seq, dx_cpt_rep_yrs, cpt_genie_sample_id)
  } else {
    rtn %<>% select(record_id, ca_seq, dx_cpt_rep_yrs, cpt_genie_sample_id)
  }
  
  return(rtn)
}
