# For features common to main GENIE and BPC - compare.

library(purrr)
library(fs)
library(here)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

char_bpc <- readr::read_rds(
  here('data', 'cohort', 'formatted_characteristics.rds')
)

# We will add in the sample year because that's comparable to main GENIE.
ca_ind <- readr::read_rds(here('data', 'cohort', 'ca_ind.rds'))
cpt <- readr::read_rds(
  here('data', 'cohort', 'cpt.rds')
)

age_at_seq_bpc <- cpt %>%
  filter(
    cpt_genie_sample_id %in%
      get_first_cpt(ca_ind, cpt, include_sample_id = T)$cpt_genie_sample_id
  ) %>%
  select(record_id, age_at_seq = dob_cpt_report_yrs)

# Double check that this is unique on record_id alone:
chk_uni_seq_age <- age_at_seq_bpc %>%
  count(record_id) %>%
  pull(n) %>%
  max %>%
  is_less_than(., 2)
if (!chk_uni_seq_age) {
  cli_abort("Non-unique sequence ages pulled from BPC.")
}

char_bpc %<>%
  left_join(., age_at_seq_bpc, by = 'record_id')

# Pull in the appropriate main GENIE data:
mg_pt <- readr::read_tsv(
  here('data-raw', 'genomic', 'main_genie', 'data_clinical_patient.txt'),
  comment = "#"
)
mg_samp <- readr::read_tsv(
  here('data-raw', 'genomic', 'main_genie', 'data_clinical_sample.txt'),
  comment = "#"
)
blad_samp_mg <- readr::read_rds(
  here('data', 'genomic', 'main_genie', 'bladder_samples_mg.rds')
)

mg_samp %<>%
  filter(SAMPLE_ID %in% blad_samp_mg) %>%
  rename_all(tolower) %>%
  select(
    record_id = patient_id,
    sample_id,
    age_at_seq_report_days,
    cancer_type_detailed
  )
mg_pt %<>%
  filter(PATIENT_ID %in% mg_samp$record_id) %>%
  rename_all(tolower) %>%
  select(
    record_id = patient_id,
    sex,
    race = primary_race,
    ethnicity,
    center
  )

# Take only the first sample for each person.
mg_samp %<>%
  mutate(
    # doing two things here - yrs to days and fixing the weird string storage system in main GENIE.
    age_at_seq = case_when(
      age_at_seq_report_days %in% ">32485" ~ 90,
      age_at_seq_report_days %in% "<6570" ~ 17,
      age_at_seq_report_days %in% "Unknown" ~ NA_real_,
      T ~ suppressWarnings(as.numeric(age_at_seq_report_days) / 365.25)
    )
  ) %>%
  select(-age_at_seq_report_days) %>%
  group_by(record_id) %>%
  arrange(age_at_seq) %>%
  slice(1) %>%
  ungroup(.)

readr::write_rds(
  mg_samp,
  here('data', 'genomic', 'main_genie', 'bladder_samples_mg_first.rds')
)

char_mg <- left_join(
  mg_pt,
  mg_samp,
  by = 'record_id',
  relationship = 'one-to-one'
)

char_mg %<>%
  mutate(
    `Sex at birth` = factor(sex),
    # Of course all the levels are different:
    `Ethnicity` = case_when(
      ethnicity %in% "Non-Spanish/non-Hispanic" ~ "Non-Spanish; non-Hispanic",
      ethnicity %in% "Spanish/Hispanic" ~ "Hispanic or Spanish",
      T ~ "Unknown whether Spanish or not"
    ),
    `Race (primary)` = case_when(
      race %in% "White" ~ "White",
      race %in% "Black" ~ "Black",
      race %in% "Asian" ~ "Asian",
      T ~ "Other or Unk."
    ),
    Institution = center,
    `Cancer type` = case_when(
      cancer_type_detailed %in% "Upper Tract Urothelial Carcinoma" ~
        "Upper tract",
      T ~ "Bladder Cancer"
    )
  ) %>%
  select(
    record_id,
    age_at_seq,
    `Sex at birth`:`Cancer type`
  )

char_comb <- bind_rows(
  mutate(study = 'BPC', char_bpc),
  mutate(study = 'Main GENIE', char_mg)
) %>%
  rename(`Age at seq.` = age_at_seq) %>%
  select(
    study,
    `Age at seq.`,
    `Sex at birth`,
    `Ethnicity`,
    `Race (primary)`,
    `Cancer type`,
    Institution,
    everything()
  )

readr::write_rds(
  char_comb,
  here('data', 'cohort', 'formatted_characteristics_with_main_genie.rds')
)
