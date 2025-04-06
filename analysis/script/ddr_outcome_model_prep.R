library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

dir_out <- here('data', 'genomic', 'ddr_def_compare')

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
  # it's weird to merge in the ca_seq and NOT the record_id from cpt_sub, since those go hand in hand.
  # instead I'll just convert the one we have, check that it's as expected, then toss it.
  rename(ddr_flags, record_id_check = record_id),
  cpt_sub,
  by = c(sample_id = 'cpt_genie_sample_id')
)
if (ddr_outcome %>% mutate(chk = record_id_check == record_id) %>% pull(chk) %>% all) {
  ddr_outcome %<>% select(-record_id_check)
} else {
  cli_abort("record_id in the ddr_outcome data does not match cpt pull-in.  check and fix.")
}

# There was one person with no path procedure date at the time this comment was written.
# I'm just going to exclude them for now.
ddr_outcome %<>%
  filter(!is.na(dx_path_proc_cpt_days))

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

rc_vec <- c('record_id', 'ca_seq')

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
  cli_abort("A status was not found for some people - check logic and fix (everyone should at least be non-invasive after dx time")
}
if (1 < (met_invasive_status %>% count(record_id, ca_seq) %>% pull(n) %>% max)) {
  cli_abort("Multiple statuses pulled - have a look at the logic/tolerances above")
}

ddr_outcome <- met_invasive_status %>%
  select(record_id, ca_seq, met_inv_status = status) %>%
  left_join(ddr_outcome, ., by = c('record_id', 'ca_seq'))

readr::write_rds(
  ddr_outcome,
  here('data', 'genomic', 'ddr_def_compare', 'ddr_outcome_premodel.rds')
)

ddr_outcome_mod <- ddr_outcome %>%
  select(sample_id, ddr, sample_type,
         dx_path_proc_cpt_days,
         de_novo_met,
         dob_ca_dx_yrs,
         upper_tract,
         institution,
         birth_year,
         race_eth,
         female,
         md_ecog_imp_num,
         met_inv_status)

ddr_outcome_mod <- ddr_outcome_mod %>% 
  rename(race = race_eth) %>% # I'm not sure why this says race_eth, just looks like race to me.
  mutate(
    met_sample = sample_type %in% "Metastasis site unspecified",
    dx_path_proc_cpt_yrs = dx_path_proc_cpt_days / 365.25, # same units.
    # This is probably best handled in an ordinal model - but it's just too much with the multiple imputation and I think it's important to get performance status in there.
    mibc_at_cpt = met_inv_status %in% c("Invasive", "Metastatic"),
    met_at_cpt = met_inv_status %in% c("Metastatic")
  ) %>%
  select(-c(sample_type,
            dx_path_proc_cpt_days,
            met_inv_status)) %>%
  replace_na(list(de_novo_met = F)) %>%
  fastDummies::dummy_cols(
    select_columns = c('institution', 'race'),
    remove_most_frequent_dummy = T,
    remove_selected_columns = T
  )

readr::write_rds(
  ddr_outcome,
  here('data', 'genomic', 'ddr_def_compare', 'ddr_outcome_mod_ready.rds')
)
    

