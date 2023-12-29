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

# Update:  This table makes little to no sense for people who never reached
#  metastatic disease.  We filter to avoid confusion.
dft_reg_neo %<>%
  filter(!is.na(dmet_time_yrs))

dft_reg_neo %<>%
  filter(!str_detect(regimen_drugs, "Investigational Drug")) %>%
  mutate(
    tt_met_reg_start_yrs = dmet_time_yrs - dx_reg_start_int_yrs,
    tt_met_reg_end_yrs = dmet_time_yrs - dx_reg_end_all_int_yrs,
    is_neoadjuvant = case_when(
      is.na(dmet_time_yrs) ~ NA,
      # this can still be NA for regimens that appeared only after met
      #  (these regimens were never presented to research team for review)
      is.na(valid_neoadjuvant) ~ NA,
      # small tolerance here - half a day - for rounding errors.
      valid_neoadjuvant %in% T & tt_met_reg_start_yrs >= -(0.5/365.25) ~ T,
      T ~ F
    )
  )

dft_multi_neo <- dft_reg_neo %>%
  filter(is_neoadjuvant) %>%
  count(record_id, ca_seq, sort = T) %>%
  filter(n > 1)
if (nrow(dft_multi_neo)) {
  cli::cli_alert_info(
    "Cases with multiple neoadjuvant regimens:"
  )
  print(dft_multi_neo)
}





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
      # judgement call here
      # if the end date is unspecified we assume it could have been within 6 months of met.
      is.na(is_gt_6mo) ~ lev_met_neoadj[1],
      is_neoadjuvant & is_gt_6mo ~ lev_met_neoadj[2],
      is_neoadjuvant & is_lte_6mo ~ lev_met_neoadj[1],
      T ~ NA_character_ # hopefully none.
    ),
    met_neoadj_f = factor(met_neoadj_f, levels = lev_met_neoadj)
  )
chk_met_neoadj_f <- dft_reg_neo %>% 
  filter(is.na(met_neoadj_f)) %>%
  nrow(.)

if (chk_met_neoadj_f > 0) {
  cli::cli_abort("There are regimens that were not categorized for 'met_neoadj_f' column.")
}



readr::write_rds(
  dft_reg_neo,
  file = here('data', 'dmet', 'neoadjuvant_status_reg.rds')
)


lev_met_neoadj_cases <- c(
  "Neoadjuvant, some ended <=6 months prior to met",
  "Neoadjuvant, all ended >6 months prior to met",
  "No neoadjuvant treatment"
)

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

# One more thing needs to be done here to make this more sensible:  We need to 
#   add back in the people who had zero regimens before met.  Then they will be 
#   correctly added as having no neoadjuvant therapy. 

dft_no_reg_addins <- dft_dmet_timing %>%
  filter(!(record_id %in% dft_reg_neo_cases$record_id)) %>%
  select(record_id, ca_seq) %>%
  mutate(
    met_neoadj_case_f = lev_met_neoadj_cases[3]
  )

dft_reg_neo_cases %<>%
  bind_rows(., dft_no_reg_addins) %>%
  mutate(
    met_neoadj_case_f = factor(met_neoadj_case_f, levels = lev_met_neoadj_cases)
  )

readr::write_rds(
  dft_reg_neo_cases,
  file = here('data', 'dmet', 'neoadjuvant_status_case.rds')
)
      

