# Description:  Creates the dataframe describing when participants had 
#   neoadjuvant treatment (if ever).

library(fs); library(purrr); library(here);
purrr::walk(.x = fs::dir_ls('R'), .f = source)


dft_ca_ind <- readr::read_rds(
  here('data', 'cohort', 'ca_ind.rds')
)
dft_reg <- readr::read_rds(
  here('data', 'cohort', 'reg.rds')
)

# This file was created by going through the list of regimens at a meeting 
#   one day.  It's totally preference of the team working on this paper.
dft_neoadj <- readr::read_csv(
  here('data', 'drug_mapping', 'neoadjuvant_regimens_2023_11_06.csv')
)

# Ensure all entries are 'yes' or 'no' only.
chk_neoadj_col <- dft_neoadj %>%
  pull(valid_neoadjuvant) %>%
  `%in%`(., c("yes", "no")) %>%
  all(.)

if (!chk_neoadj_col) {
  cli::cli_abort("Invalid entries in the input neoadjuvant CSV (valid_neoadjuvant column)")
}

# Make sure all regimens are in the data - i.e. no typos or corrupt regimens.
chk_reg_drug_col <- dft_neoadj %>%
  pull(regimen_drugs) %>%
  `%in%`(
    .,
    unique(pull(dft_reg, regimen_drugs))
  ) %>%
  all(.)

if (!chk_reg_drug_col) {
  cli::cli_abort("Invalid entries in the the input neoadjuvant CSV (regimen_drugs column)")
}

dft_neoadj %<>%
  arrange(regimen_drugs) %>%
  mutate(
    valid_neoadjuvant = case_when(
      valid_neoadjuvant %in% "yes" ~ T,
      valid_neoadjuvant %in% "no" ~ F,
      T ~ NA
    )
  )

# write a copy of this:
readr::write_rds(
  x = dft_neoadj,
  file = here('data', 'drug_mapping', 'neoadjuvant_reg_key.rds')
)


dft_reg_neo <- dft_reg %>% 
  left_join(
    ., 
    dft_neoadj,
    by = "regimen_drugs"
  )

dft_reg_neo %<>%
  select(
    record_id,
    ca_seq,
    regimen_number,
    regimen_drugs,
    valid_neoadjuvant,
    dx_reg_start_int_yrs,
    dx_reg_end_all_int_yrs
  )

dft_dmet_timing <- make_dmet_status_block(dft_ca_ind) %>%
  filter(dmet_status %in% "Distant Metastasis") %>%
  select(record_id, ca_seq, dmet_time_yrs = dx_block_start)

dft_reg_neo <- left_join(
  dft_reg_neo,
  dft_dmet_timing,
  by = c("record_id", "ca_seq")
)

dft_reg_neo %<>%
  filter(!str_detect(regimen_drugs, "Investigational Drug")) %>%
  mutate(
    tt_met_reg_start_yrs = dmet_time_yrs - dx_reg_start_int_yrs,
    tt_met_reg_end_yrs = dmet_time_yrs - dx_reg_end_all_int_yrs,
    is_neoadjuvant = case_when(
      is.na(dmet_time_yrs) ~ NA,
      # this can still be NA for regimens that appeared only after met
      #. (these regimens were never presented to physicians for review)
      is.na(valid_neoadjuvant) ~ NA,
      # small tolerance here - half a day - for rounding errors.
      valid_neoadjuvant %in% T & tt_met_reg_start_yrs >= -(0.5/365.25) ~ T,
      T ~ F
    )
  )

lev_met_neoadj <- c(
  "Neoadjuvant, ended <=6 months prior to met",
  "Neoadjuvant, ended >6 months prior to met",
  "Not neoadjuvant treatment"
)

dft_reg_neo %<>%
  mutate(
    is_gt_6mo = tt_met_reg_end_yrs > 0.5,
    is_lte_6mo = tt_met_reg_end_yrs <= 0.5,
    met_neoadj_f = case_when(
      is.na(is_neoadjuvant) ~ lev_met_neoadj[3],
      !is_neoadjuvant ~ lev_met_neoadj[3],
      is_neoadjuvant & is_gt_6mo ~ lev_met_neoadj[2],
      is_neoadjuvant & is_lte_6mo ~ lev_met_neoadj[1],
      # These are backup designations. If we don't have the end date, we go
      #   with what the start date tells us.
      is_neoadjuvant & (tt_met_reg_start_yrs <= 0.5) ~ lev_met_neoadj[1],
      T ~ lev_met_neoadj[2], # definitely debatable.
    ),
    met_neoadj_f = factor(met_neoadj_f, levels = lev_met_neoadj)
  )

readr::write_rds(
  dft_reg_neo,
  file = here('data', 'dmet', 'neoadjuvant_status_reg.rds')
)

# The case labels are similar, just need to update the last one:
lev_met_neoadj_cases <- lev_met_neoadj
lev_met_neoadj_cases[3] <- "No neoadjuvant treatment"

dft_reg_neo_cases <- dft_reg_neo %<>%
  group_by(record_id, ca_seq) %>%
  summarize(
    met_neoadj_case_f = case_when(
      any(met_neoadj_f %in% lev_met_neoadj[1]) ~ lev_met_neoadj_cases[1],
      any(met_neoadj_f %in% lev_met_neoadj[2]) ~ lev_met_neoadj_cases[2],
      T ~ lev_met_neoadj_cases[3]
    ),
    .groups = "drop"
  )

readr::write_rds(
  dft_reg_neo_cases,
  file = here('data', 'dmet', 'neoadjuvant_status_case.rds')
)
      

# Todo:
#  - update the manuscript sections to have 3 tables.
#  - remove repeats in the tables if not already done.
#  - remove same-class regimens in the table (e.g. carboplatin and cisplatin)
