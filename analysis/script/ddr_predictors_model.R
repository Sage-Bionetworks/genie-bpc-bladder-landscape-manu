library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

# These DDR flags only cover each person's first sample in BPC:
ddr_flags <- readr::read_rds(
  here('data', 'genomic', 'ddr_def_compare', 'ddr_flags_first_sample.rds')
)
cpt <- readr::read_rds(here('data', 'cohort', "cpt_aug.rds"))
timeless_features <- readr::read_rds(
  here('data', 'cohort', 'time_invariant_model_factors.rds')
)
med_onc <- readr::read_rds(here('data', 'cohort', "med_onc.rds"))
med_onc %<>% augment_med_onc_imputed_ecog(., add_numeric = T)



cpt_sub <- cpt %>% select(
  cpt_genie_sample_id, record_id, ca_seq, 
  sample_type,
  dx_path_proc_cpt_days,
  dob_cpt_report_days, dx_cpt_rep_days
) %>%
  mutate(
    dob_dx_days = dob_cpt_report_days - dx_cpt_rep_days,
    dob_path_proc_cpt_days = dx_path_proc_cpt_days + dob_dx_days
  ) %>%
  select(-c(dob_cpt_report_days, dx_cpt_rep_days))

# ddr_outcome = DDR as an outcome variable - distinguishing from most other uses we have.
ddr_outcome <- left_join(
  ddr_flags,
  cpt_sub,
  by = c(sample_id = 'cpt_genie_sample_id')
) 

