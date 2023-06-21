get_first_regimen <- function(ca_ind_dat, reg_dat) {
  rtn <- ca_ind_dat %>% 
    select(record_id, ca_seq) %>%
    left_join(., reg_dat, by = c("record_id", "ca_seq")) %>%
    arrange(dx_reg_start_int_yrs) %>%
    group_by(record_id) %>%
    slice(1) %>%
    ungroup() %>%
    select(record_id, ca_seq, dx_reg_start_int_yrs)
  
  return(rtn)
}
