# This is a similar script to derive data for DDR analyses of 1L therapy.
# The main difference is we want ALL 1L therapies here, with flags for DDR
#   and flags for platinum.

library(purrr)
library(here)
library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source) # also load

dft_pt <- readr::read_rds(here('data', 'cohort', "pt.rds"))
dft_ca_ind <- readr::read_rds(here('data', 'cohort', "ca_ind.rds"))
dft_reg <- readr::read_rds(here('data', 'cohort', "reg.rds"))
dft_cpt <- readr::read_rds(here('data', 'cohort', "cpt_aug.rds"))
dft_alt <- readr::read_rds(here('data', 'genomic', 'alterations.rds'))
dft_med_onc <- readr::read_rds(here('data', 'cohort', "med_onc.rds"))
dft_med_onc %<>%
  augment_med_onc_imputed_ecog(., add_numeric = T)
lot <- readr::read_rds(
  here('data', 'dmet', 'lines_of_therapy', 'lot.rds')
)

dir_out <- here('data', 'survival', 'ddr_onco_1L')


rel_reg <- lot %>%
  filter(line_therapy %in% 1) %>%
  left_join(
    .,
    select(
      dft_reg,
      record_id,
      ca_seq,
      contains('regimen_number'),
      dx_reg_start_int,
      dx_reg_end_all_int,
      os_g_status,
      tt_os_g_days
    ),
    by = c('record_id', 'ca_seq', 'regimen_number'),
    relationship = 'one-to-one'
  )

# Add in the regimen start time from dob.
# Have to reconstruct as it's not available in this cohort.
rec_ca <- c('record_id', 'ca_seq')
rel_reg <- left_join(
  rel_reg,
  select(dft_ca_ind, record_id, ca_seq, dob_ca_dx_days),
  by = rec_ca
) %>%
  mutate(
    dob_reg_start_days = dx_reg_start_int + dob_ca_dx_days
  ) %>%
  relocate(dob_reg_start_days, .before = dx_reg_start_int)


rel_reg <- add_entry_time(
  dat_index = rel_reg,
  var_index = "dx_reg_start_int",
  dft_ca_ind,
  dft_cpt
)


# Add in the genomic covariates:
custom_ddr_list <- c(
  'ERCC2',
  'ERCC5',
  'BRCA1',
  'BRCA2',
  'RECQL4',
  'RAD51C',
  'ATM',
  'ATR',
  'FANCC'
) %>%
  sort(.)
dft_onco_ddr <- dft_alt %>%
  filter(hugo %in% custom_ddr_list) %>%
  filter(oncogenic %in% c("Likely Oncogenic", "Oncogenic"))
# Add in the relevant stuff from CPT data:
dft_onco_ddr <- dft_cpt %>%
  select(
    record_id,
    ca_seq,
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
  left_join(
    .,
    select(rel_reg, record_id, ca_seq, dx_entry),
    by = c('record_id', 'ca_seq')
  ) %>%
  # take only those who had a path proc time before entry.
  filter(dx_entry > dx_path_proc_cpt_days - 0.5)
rel_reg <- onco_ddr_before_entry %>%
  group_by(record_id, ca_seq) %>%
  summarize(ddr_before_entry = T) %>%
  left_join(rel_reg, ., by = c('record_id', 'ca_seq')) %>%
  replace_na(list(ddr_before_entry = F))


rel_reg <- rel_reg %>%
  mutate(
    regimen_cat = categorize_lines(regimen_drugs)
  ) %>%
  relocate(regimen_cat, .after = regimen_drugs)


# Adding a few things here which might be valuable as confounders:
dft_extra_var <- readr::read_rds(
  here('data', 'cohort', 'time_invariant_model_factors.rds')
) %>%
  filter(!is.na(dx_dmet_yrs)) %>%
  select(
    record_id,
    ca_seq,
    dx_dmet_yrs,
    institution,
    birth_year,
    race_eth,
    female
  )
rel_reg <- rel_reg %>%
  left_join(., dft_extra_var, by = c('record_id', 'ca_seq'))

# Add in the latest med onc observation:
dft_latest_med_onc <- get_med_onc_by_timing(
  dat_med_onc = dft_med_onc,
  dat_cutoff = select(rel_reg, record_id, dob_reg_start_days),
  var_cutoff = "dob_reg_start_days",
  remove_missing = T
)
rel_reg %<>%
  left_join(
    .,
    select(dft_latest_med_onc, record_id, md_ecog_imp_num),
    by = "record_id"
  )

# Add in de novo met indicator:
rel_reg <- get_dmet_time(dft_ca_ind, annotate_type = T) %>%
  mutate(
    de_novo_met = case_when(
      .met_type %in% "stage_iv_with_immediate_dmet" ~ T,
      T ~ F
    )
  ) %>%
  select(record_id, ca_seq, de_novo_met) %>%
  left_join(
    rel_reg,
    .,
    by = c('record_id', 'ca_seq')
  )

rel_reg <- rel_reg %>%
  # int is less clear than days
  rename(
    dx_reg_start_days = dx_reg_start_int,
    dx_reg_end_all_days = dx_reg_end_all_int,
    dx_entry_days = dx_entry
  ) %>%
  mutate(reg_fcpt_days = dx_first_cpt_rep_days - dx_reg_start_days)

# Convert to years and locate the similar variables together:
rel_reg <- rel_reg %>%
  # obviously a great opportunity for a "covert and rename all these vars" function...
  mutate(
    dob_reg_start_yrs = day_to_year_genie(dob_reg_start_days),
    dx_reg_start_yrs = day_to_year_genie(dx_reg_start_days),
    dx_reg_end_all_yrs = day_to_year_genie(dx_reg_end_all_days),
    tt_os_g_yrs = day_to_year_genie(tt_os_g_days),
    dob_ca_dx_yrs = day_to_year_genie(dob_ca_dx_days),
    dx_first_cpt_rep_yrs = day_to_year_genie(dx_first_cpt_rep_days),
    dx_entry_yrs = day_to_year_genie(dx_entry_days),
    reg_fcpt_yrs = day_to_year_genie(reg_fcpt_days)
  ) %>%
  relocate(dob_reg_start_yrs, .after = dob_reg_start_days) %>%
  relocate(dx_reg_start_yrs, .after = dx_reg_start_days) %>%
  relocate(dx_reg_end_all_yrs, .after = dx_reg_end_all_days) %>%
  relocate(tt_os_g_yrs, .after = tt_os_g_days) %>%
  relocate(dob_ca_dx_yrs, .after = dob_ca_dx_days) %>%
  relocate(dx_first_cpt_rep_yrs, .after = dx_first_cpt_rep_days) %>%
  relocate(dx_entry_yrs, .after = dx_entry_days) %>%
  relocate(reg_fcpt_yrs, .after = reg_fcpt_days)


readr::write_rds(
  rel_reg,
  here(dir_out, 'met_ddr_surv_all_1L.rds')
)


# Then to make the cisplatin/carboplatin only version is easy:
plat_only <- rel_reg %>%
  filter(str_detect(regimen_drugs, "Carboplat|Cisplat"))

readr::write_rds(
  plat_only,
  here(dir_out, 'met_ddr_surv_plat_1L.rds')
)
