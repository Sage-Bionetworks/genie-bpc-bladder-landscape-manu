# This is a derivation of time-invariant factors which may be useful in modeling.  The idea is to save time on data cleaning steps which will be 
# the same over many models.
# I will just keep expanding this each time it makes sense, so in practice
#   it would be best to load, subset, merge this (subset being key).

library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

dft_pt <- readr::read_rds(here('data', 'cohort', "pt.rds"))
dft_ca_ind <- readr::read_rds(here('data', 'cohort', "ca_ind.rds"))
dft_reg <- readr::read_rds(here('data', 'cohort', "reg.rds"))

dft_met_timing <- get_dmet_time(dft_ca_ind)

rc_vec <- c('record_id', 'ca_seq') # space saver.
dft_time_invar <- dft_ca_ind %>%
  select(record_id, ca_seq)

dft_time_invar <- dft_met_timing %>%
  mutate(de_novo_met = if_else(abs(dx_dmet_yrs) < 0.0005, T, F)) %>%
  select(record_id, ca_seq, dx_dmet_yrs, de_novo_met) %>%
  left_join(dft_time_invar, ., by = rc_vec)

dft_time_invar <- dft_ca_ind %>% 
  select(record_id, ca_seq, dob_ca_dx_yrs) %>%
  left_join(dft_time_invar, ., by = rc_vec) 

dft_time_invar %<>%
  left_join(
    .,
    select(dft_pt, record_id, institution, birth_year, naaccr_race_code_primary, naaccr_ethnicity_code, naaccr_sex_code),
    by = "record_id"
  ) 

dft_time_invar %<>%
  mutate(
    race_eth = case_when(
      # don't think we can support this AND a term for UHN/MSK.
      # naaccr_ethnicity_code %in% c("Puerto Rican", "Spanish NOS or Hispanic NOS or Latino NOS") ~ "Hispanic",
      naaccr_race_code_primary %in% c("White") ~ "White",
      naaccr_race_code_primary %in% c("Black") ~ "Black",
      
      naaccr_race_code_primary %in% c("Chinese", "Other Asian", "Asian Indian or Pakistani NOS") ~ "Asian",
      T ~ "race_other_unk"
    ),
    female = if_else(naaccr_sex_code %in% "Female", T, F)
  ) %>%
  select(-naaccr_sex_code, -naaccr_ethnicity_code, -naaccr_race_code_primary)

readr::write_rds(
  x = dft_time_invar,
  file = here('data', 'cohort', 'time_invariant_model_factors.rds')
)
  


# 
# 
# 
# dft_reg_start_dob <- left_join(
#   select(dft_met_ddr_surv, record_id, ca_seq, dx_reg_start_int_yrs),
#   select(dft_ca_ind, record_id, ca_seq, dob_ca_dx_yrs),
#   by = rc_vec
# ) %>%
#   mutate(
#     dob_reg_start_yrs = dx_reg_start_int_yrs + dob_ca_dx_yrs,
#     dob_reg_start_days = dob_reg_start_yrs * 365.25
#   ) %>%
#   select(record_id, dob_reg_start_days) 
# 
# dft_time_invar %<>%
#   left_join(., dft_reg_start_dob, by = 'record_id')
# 
# dft_latest_med_onc <- get_med_onc_by_timing(
#   dat_med_onc = dft_med_onc,
#   dat_cutoff = dft_reg_start_dob,
#   var_cutoff = "dob_reg_start_days",
#   remove_missing = T
# )
# 
# dft_extra_var %<>%
#   left_join(
#     .,
#     select(dft_latest_med_onc, record_id, md_ecog_imp_num),
#     by = "record_id"
#   )
