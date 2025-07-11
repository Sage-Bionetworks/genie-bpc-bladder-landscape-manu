library(purrr)
library(here)
library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

dir_out <- here('data', 'genomic', 'ddr_def_compare', 'ddr_as_outcome')

# These DDR flags only cover each person's first sample in BPC:
ddr_flags <- readr::read_rds(
  here('data', 'genomic', 'ddr_def_compare', 'ddr_flags_first_sample.rds')
)
ddr_flags %<>% rename(record_id = patient_id)
ca_ind <- readr::read_rds(here('data', 'cohort', 'ca_ind.rds'))
cpt <- readr::read_rds(here('data', 'cohort', "cpt_aug.rds"))
timeless_features <- readr::read_rds(
  here('data', 'cohort', 'time_invariant_model_factors.rds')
)
med_onc <- readr::read_rds(here('data', 'cohort', "med_onc.rds"))
med_onc %<>% augment_med_onc_imputed_ecog(., add_numeric = T)


cpt_sub <- cpt %>%
  select(
    cpt_genie_sample_id,
    record_id,
    ca_seq,
    sample_type,
    dx_path_proc_cpt_days,
    dob_cpt_report_days,
    dx_cpt_rep_days
  ) %>%
  mutate(
    dob_dx_days = dob_cpt_report_days - dx_cpt_rep_days,
    dob_path_proc_cpt_days = dx_path_proc_cpt_days + dob_dx_days
  ) %>%
  select(-c(dob_cpt_report_days, dx_cpt_rep_days))

# ddr_outcome = DDR as an outcome variable - distinguishing from most other uses we have.
ddr_outcome <- left_join(
  # it's weird to merge in the ca_seq and NOT the record_id from cpt_sub, since those go hand in hand.
  # instead I'll just convert the one we have, check that it's as expected, then toss it.
  rename(ddr_flags, record_id_check = record_id),
  cpt_sub,
  by = c(sample_id = 'cpt_genie_sample_id')
)
if (
  ddr_outcome %>%
    mutate(chk = record_id_check == record_id) %>%
    pull(chk) %>%
    all
) {
  ddr_outcome %<>% select(-record_id_check)
} else {
  cli_abort(
    "record_id in the ddr_outcome data does not match cpt pull-in.  check and fix."
  )
}

# There was one person with no path procedure date at the time this comment was written.
# I'm just going to exclude them for now.
ddr_outcome %<>%
  filter(!is.na(dx_path_proc_cpt_days))

rc_vec <- c('record_id', 'ca_seq')
ddr_outcome <- left_join(
  ddr_outcome,
  timeless_features,
  by = rc_vec
)


# Get time varying model factors - pulling the last knowledge before CPT procedure.
latest_med_onc <- get_med_onc_by_timing(
  dat_med_onc = med_onc,
  dat_cutoff = select(ddr_outcome, record_id, dob_path_proc_cpt_days),
  var_cutoff = "dob_path_proc_cpt_days",
  remove_missing = T
)
ddr_outcome %<>%
  left_join(
    .,
    select(latest_med_onc, record_id, md_ecog_imp_num),
    by = "record_id"
  )

met_invasive_status <- make_dmet_musc_prop_status_block(ca_ind) %>%
  mutate(
    dx_block_start_days = dx_block_start * 365.25,
    dx_block_end_days = dx_block_end * 365.25
  ) %>%
  select(-c(dx_block_start, dx_block_end)) %>%
  left_join(
    .,
    select(ddr_outcome, record_id, ca_seq, dx_path_proc_cpt_days),
    by = rc_vec
  )

# filter down to the status only at the time of sampling.
met_invasive_status <- met_invasive_status %>%
  filter(
    dx_block_start_days <= dx_path_proc_cpt_days + 0.5, # half day tolerance for rounding.
    dx_block_end_days >= dx_path_proc_cpt_days + 0.5 # if it's on the border we take the newly declared one on the procedure date.
  )

if (nrow(anti_join(ddr_outcome, met_invasive_status, by = rc_vec)) > 0) {
  cli_abort(
    "A status was not found for some people - check logic and fix (everyone should at least be non-invasive after dx time"
  )
}
if (
  1 < (met_invasive_status %>% count(record_id, ca_seq) %>% pull(n) %>% max)
) {
  cli_abort(
    "Multiple statuses pulled - have a look at the logic/tolerances above"
  )
}

ddr_outcome <- met_invasive_status %>%
  select(record_id, ca_seq, met_inv_status = status) %>%
  left_join(ddr_outcome, ., by = c('record_id', 'ca_seq'))

panel_size_by_sample <- readr::read_rds(
  here('data', 'genomic', 'main_genie', 'panel_size_by_sample.rds')
)

ddr_outcome <- left_join(
  ddr_outcome,
  select(panel_size_by_sample, sample_id, seq_assay_id, n_genes),
  by = 'sample_id'
)

readr::write_rds(
  ddr_outcome,
  here(dir_out, 'ddr_outcome_premodel.rds')
)

ddr_outcome_mod <- ddr_outcome %>%
  select(
    sample_id,
    ddr,
    sample_type,
    dx_path_proc_cpt_days,
    de_novo_met,
    dob_ca_dx_yrs,
    upper_tract,
    institution,
    birth_year,
    race_eth,
    female,
    md_ecog_imp_num,
    met_inv_status,
    seq_assay_id,
    n_genes
  )

ddr_outcome_mod <- ddr_outcome_mod %>%
  rename(race = race_eth) %>% # I'm not sure why this says race_eth, just looks like race to me.
  mutate(
    met_sample = sample_type %in% "Metastasis site unspecified",
    dx_path_proc_cpt_yrs = dx_path_proc_cpt_days / 365.25, # same units.
    # This is probably best handled in an ordinal model - but it's just too much with the multiple imputation and I think it's important to get performance status in there.
    mibc_at_cpt = met_inv_status %in% c("Invasive", "Metastatic"),
    met_at_cpt = met_inv_status %in% c("Metastatic")
  ) %>%
  select(-c(sample_type, dx_path_proc_cpt_days, met_inv_status)) %>%
  replace_na(list(de_novo_met = F)) %>%
  fastDummies::dummy_cols(
    select_columns = c('institution', 'race'),
    remove_most_frequent_dummy = T,
    remove_selected_columns = T
  )

ddr_outcome_mod %<>%
  mutate(
    panel_genes_100 = n_genes / 100,
    birth_year_10 = birth_year / 10,
    age_dx_10 = dob_ca_dx_yrs / 10
  ) %>%
  select(-c(n_genes, birth_year, dob_ca_dx_yrs))


readr::write_rds(
  ddr_outcome_mod,
  here(dir_out, 'ddr_outcome_mod_ready.rds')
)


# Finally try to do the same for main GENIE - or as close as we can get.
ddr_flags_main <- readr::read_rds(
  here(
    'data',
    'genomic',
    'ddr_def_compare',
    'ddr_flags_first_sample_main_genie.rds'
  )
)

sample_main <- readr::read_tsv(
  here('data-raw', 'genomic', 'main_genie', 'data_clinical_sample.txt'),
  skip = 4
) %>%
  rename_all(tolower)

ddr_flags_main <- sample_main %>%
  select(
    sample_id,
    age_at_seq_report_days,
    oncotree_code,
    sample_type,
    cancer_type,
    seq_year
  ) %>%
  left_join(
    ddr_flags_main,
    .,
    by = 'sample_id'
  )
pt_main <- readr::read_tsv(
  here('data-raw', 'genomic', 'main_genie', 'data_clinical_patient.txt'),
  skip = 4
) %>%
  rename_all(tolower)
ddr_flags_main <- pt_main %>%
  select(
    patient_id,
    sex,
    primary_race,
    birth_year,
    center
  ) %>%
  left_join(
    ddr_flags_main,
    .,
    by = 'patient_id'
  )

ddr_flags_main <- left_join(
  ddr_flags_main,
  select(panel_size_by_sample, sample_id, seq_assay_id, n_genes),
  by = 'sample_id'
)

ddr_flags_main %<>%
  mutate(
    age_seq_yrs = case_when(
      is.na(age_at_seq_report_days) | age_at_seq_report_days %in% "Unknown" ~
        NA_real_,
      age_at_seq_report_days %in% ">32485" ~ 90, # could be higher - this is fine.
      age_at_seq_report_days %in% "<6570" ~ 17, # could be lower - this is fine.
      T ~ as.numeric(age_at_seq_report_days) / 365.25
    ),
    birth_year = case_when(
      birth_year %in% "Unknown" ~ NA_real_,
      birth_year %in% "cannotReleaseHIPAA" ~ 1925, # could be higher, this is OK.
      birth_year %in% "withheld" ~ 2005, # could be lower, this is OK.
      T ~ as.numeric(birth_year)
    )
  ) %>%
  select(-age_at_seq_report_days)

# K there's some obvious problems here:
# ggplot(ddr_flags_main, aes(x = age_seq_yrs, y = birth_year)) + geom_point()

ddr_flags_main %<>%
  filter(
    !(age_seq_yrs > 40 & birth_year > 2000)
  )

# ggplot(ddr_flags_main, aes(x = age_seq_yrs, y = birth_year)) + geom_point()
# Good enough - probably still some impossible things but at least the mega high leverage points are gone.

ddr_flags_main %<>%
  mutate(
    upper_tract = oncotree_code %in% "UTUC",
    met_sample = sample_type %in% "Metastasis",
    female = sex %in% "Female",
    race = case_when(
      primary_race %in% c("White", "Black", "Asian") ~ primary_race,
      T ~ "race_oth_unk"
    )
  ) %>%
  select(
    -c(oncotree_code, sample_type, seq_year, sex, cancer_type, primary_race)
  ) %>%
  rename(institution = center) # have seq year coded between birth_year and age at seq.

ddr_flags_main_mod <- ddr_flags_main %<>%
  fastDummies::dummy_cols(
    select_columns = c('institution', 'race'),
    remove_most_frequent_dummy = T,
    remove_selected_columns = T
  )

ddr_flags_main_mod %<>%
  mutate(
    panel_genes_100 = n_genes / 100,
    birth_year_10 = birth_year / 10,
    age_seq_10 = age_seq_yrs / 10
  ) %>%
  select(-c(n_genes, birth_year, age_seq_yrs))

readr::write_rds(
  ddr_flags_main_mod,
  here(dir_out, 'ddr_outcome_mod_ready_main.rds')
)
