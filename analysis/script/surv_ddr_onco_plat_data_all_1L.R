# This is a similar script to derive data for DDR analyses of 1L therapy.
# The main difference is we want ALL 1L therapies here, with flags for DDR
#   and flags for platinum.

library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source) # also load

dft_pt <- readr::read_rds(here('data', 'cohort', "pt.rds"))
dft_ca_ind <- readr::read_rds(here('data', 'cohort', "ca_ind.rds"))
dft_reg <- readr::read_rds(here('data', 'cohort', "reg.rds"))
dft_cpt <- readr::read_rds(here('data', 'cohort', "cpt_aug.rds"))
dft_alt <- readr::read_rds(here('data', 'genomic','alterations.rds'))
dft_med_onc <- readr::read_rds(here('data', 'cohort', "med_onc.rds"))
dft_med_onc %<>%
  augment_med_onc_imputed_ecog(., add_numeric = T)
lot <- readr::read_rds(
  here('data', 'dmet', 'lines_of_therapy', 'lot.rds')
)

dir_out <- here('data', 'survival', 'ddr_onco')




rel_reg <- lot %>%
  filter(line_therapy %in% 1) %>%
  left_join(
    .,
    select(dft_reg, record_id, ca_seq, contains('regimen_number'),
           dx_reg_start_int, dx_reg_end_all_int,
           os_g_status, tt_os_g_days),
    by = c('record_id', 'ca_seq', 'regimen_number'),
    relationship = 'one-to-one'
  )

rel_reg <- add_entry_time(
  dat_index = rel_reg,
  var_index = "dx_reg_start_int",
  dft_ca_ind,
  dft_cpt
)
  

# Add in the genomic covariates:
custom_ddr_list <- c(
  'ERCC2', 'ERCC5', 
  'BRCA1', 'BRCA2', 'RECQL4', 'RAD51C', 'ATM', 
  'ATR', 'FANCC'
) %>% sort(.)
dft_onco_ddr <- dft_alt %>% 
  filter(hugo %in% custom_ddr_list) %>%
  filter(oncogenic %in% c("Likely Oncogenic", "Oncogenic"))
# Add in the relevant stuff from CPT data:
dft_onco_ddr <- dft_cpt %>%
  select(
    record_id, ca_seq,
    dx_path_proc_cpt_days, # interval from dx to pathology procedure of this CPT. 
    dx_cpt_rep_days, # interval from dx to report date of this CPT.
    cpt_genie_sample_id
  ) %>%
  left_join(
    dft_onco_ddr,
    .,
    by = c(sample_id = "cpt_genie_sample_id")
  )
onco_ddr_before_entry <- dft_onco_ddr %>%
  left_join(., select(rel_reg, record_id, ca_seq, dx_entry),
            by = c('record_id', 'ca_seq')) %>%
  # take only those who had a path proc time before entry.
  filter(dx_entry > dx_path_proc_cpt_days - 0.5)
rel_reg <- onco_ddr_before_entry %>%
  group_by(record_id, ca_seq) %>%
  summarize(ddr_before_entry = T) %>%
  left_join(rel_reg, ., by = c('record_id', 'ca_seq')) %>%
  replace_na(list(ddr_before_entry = F)) 




cli_abort('Stopped here - I think we"re close just need a bit more cleaning')



# Two Options here.

# Option 1: To just select the first regimen, no matter how it falls in lines:
# dft_met_ddr_surv <- rel_reg %>%
#   group_by(record_id, ca_seq) %>%
#   arrange(dx_reg_start_int_yrs) %>%
#   slice(1) %>%
#   ungroup(.)

# Option 2: To select only first line therapy:
dft_lot <- readr::read_rds(here('data', 'dmet', 'lines_of_therapy', 'lot.rds'))
dft_met_ddr_surv <- dft_lot %>%
  filter(
    line_therapy %in% 1,
    str_detect(regimen_drugs, "Cisplatin|Carboplatin")
  ) %>%
  select(record_id, ca_seq, regimen_number) %>%
  left_join(
    .,
    rel_reg,
    by = c('record_id', 'ca_seq', 'regimen_number')
  )




dft_met_ddr_surv <- left_join(
  dft_met_ddr_surv,
  dft_onco_ddr_flags,
  by = c("record_id", "ca_seq")
) 

dft_met_ddr_surv %<>%
  replace_na(
    list(
      ddr_before_pm_reg = 0,
      ddr_before_entry = 0
    )
  ) 

# Still need the cohort entry time for survival to work out:
dft_met_ddr_surv %<>%
  left_join(
    ., 
    dft_first_cpt,
    by = c('record_id', 'ca_seq')
  ) %>%
  # I just cant stand these names...
  rename(
    os_first_met_reg_status = os_g_status,
    tt_os_first_met_reg_yrs = tt_os_g_yrs
  ) %>%
  mutate(
    fmr_fcpt_yrs = dx_first_cpt_rep_yrs - dx_reg_start_int_yrs
  )





# Adding a few things here which might be valuable as confounders:
dft_extra_var <- readr::read_rds(
  here('data', 'cohort', 'time_invariant_model_factors.rds')
) %>%
  filter(!is.na(dx_dmet_yrs)) %>%
  select(record_id, ca_seq, dx_dmet_yrs, dob_ca_dx_yrs, institution,
         birth_year, race_eth, female)

rc_vec <- c('record_id', 'ca_seq')
dft_reg_start_dob <- left_join(
  select(dft_met_ddr_surv, record_id, ca_seq, dx_reg_start_int_yrs),
  select(dft_ca_ind, record_id, ca_seq, dob_ca_dx_yrs),
  by = rc_vec
) %>%
  mutate(
    dob_reg_start_yrs = dx_reg_start_int_yrs + dob_ca_dx_yrs,
    dob_reg_start_days = dob_reg_start_yrs * 365.25
  ) %>%
  select(record_id, dob_reg_start_days) 

dft_extra_var %<>%
  left_join(., dft_reg_start_dob, by = 'record_id')

dft_latest_med_onc <- get_med_onc_by_timing(
  dat_med_onc = dft_med_onc,
  dat_cutoff = dft_reg_start_dob,
  var_cutoff = "dob_reg_start_days",
  remove_missing = T
)

dft_extra_var %<>%
  left_join(
    .,
    select(dft_latest_med_onc, record_id, md_ecog_imp_num),
    by = "record_id"
  )

dft_met_ddr_surv <- left_join(
  dft_met_ddr_surv,
  select(dft_extra_var, -dx_dmet_yrs), # no need to double up.
  by = rc_vec
)


readr::write_rds(
  dft_met_ddr_surv,
  here(dir_out, 'met_ddr_surv.rds')
)


















