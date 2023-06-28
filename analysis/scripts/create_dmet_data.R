# Description: Create versions of the datasets filtered down to events after metastasis.

read_wrap <- function(p) {
  read_csv(file = here("data", 'cohort', p), show_col_types = F)
}

dft_pt <- read_wrap("pt.csv")
dft_ca_ind <- read_wrap("ca_ind.csv")
dft_reg <- read_wrap("reg.csv")

dft_dmet_timing <- make_dmet_status_block(dft_ca_ind) %>%
  filter(dmet_status %in% "Distant Metastasis") %>%
  select(record_id, ca_seq, dmet_time_yrs = dx_block_start)

dft_pt_dmet <- dft_pt %>% filter(record_id %in% dft_dmet_timing$record_id)

dft_ca_ind_dmet <- left_join(
  dft_ca_ind,
  dft_dmet_timing,
  by = c("record_id", "ca_seq")
) %>%
  filter(!is.na(dmet_time_yrs)) %>%
  select(-dmet_time_yrs)

# s = "Start of regimen after distant metastasis", which we will save as below.
dft_reg_dmet_s <- left_join(
  dft_reg,
  dft_dmet_timing,
  by = c("record_id", "ca_seq")
) %>%
  filter(!is.na(dmet_time_yrs)) %>%
  filter(dx_reg_start_int_yrs >= dmet_time_yrs) %>%
  select(-dmet_time_yrs)


# Todo:  Move the drug code over here.



fs::dir_create('data', 'dmet')

write_wrap <- function(obj, file_name) {
  readr::write_csv(x = obj, file = here('data', 'dmet', file_name))
}

write_wrap(dft_pt_dmet, file_name = "pt_dmet.csv")
write_wrap(dft_ca_ind_dmet, file_name = "ca_ind_dmet.csv")
write_wrap(dft_reg_dmet_s, file_name = "reg_start_gte_dmet.csv")


