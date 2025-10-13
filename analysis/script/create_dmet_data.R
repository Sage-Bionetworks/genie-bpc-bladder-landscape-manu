# Description: Create versions of the datasets filtered down to events after metastasis.

library(purrr)
library(here)
library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

read_wrap <- function(p) {
  read_rds(
    file = here("data", 'cohort', paste0(p, ".rds"))
  )
}

dft_pt <- read_wrap("pt")
dft_ca_ind <- read_wrap("ca_ind")
dft_reg <- read_wrap("reg")


dft_dmet_timing <- get_dmet_time(dft_ca_ind)

dft_pt_dmet <- dft_pt %>% filter(record_id %in% dft_dmet_timing$record_id)

dft_ca_ind_dmet <- left_join(
  dft_ca_ind,
  dft_dmet_timing,
  by = c("record_id", "ca_seq")
) %>%
  filter(!is.na(dx_dmet_yrs)) %>%
  select(-dx_dmet_yrs)

# s = "Start of regimen after distant metastasis", which we will save as below.
dft_reg_dmet_s <- left_join(
  dft_reg,
  dft_dmet_timing,
  by = c("record_id", "ca_seq")
) %>%
  filter(!is.na(dx_dmet_yrs)) %>%
  filter(dx_reg_start_int_yrs >= dx_dmet_yrs) %>%
  select(-dx_dmet_yrs)


# Todo:  Move the drug code over here.

write_wrap <- function(obj, file_name) {
  readr::write_rds(
    x = obj,
    file = here('data', 'dmet', paste0(file_name, '.rds'))
  )
}

write_wrap(dft_pt_dmet, file_name = "pt_dmet")
write_wrap(dft_ca_ind_dmet, file_name = "ca_ind_dmet")
write_wrap(dft_reg_dmet_s, file_name = "reg_start_gte_dmet")
