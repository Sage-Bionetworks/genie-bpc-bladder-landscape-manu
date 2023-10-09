lot_filter_prep <- function(dat_reg, max_reg) {
  dat_reg %<>%
    arrange(record_id, ca_seq, regimen_number) %>%
    group_by(record_id, ca_seq) %>%
    mutate(.conseq_reg_number = 1:n()) %>%
    ungroup(.) %>%
    filter(.conseq_reg_number %in% 1:max_reg)
  return(dat_reg)
}
