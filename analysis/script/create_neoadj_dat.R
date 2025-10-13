# Description:  Creates the dataframe describing when participants had
#   neoadjuvant treatment (if ever).

library(fs)
library(purrr)
library(here)
purrr::walk(.x = fs::dir_ls('R'), .f = source)


dft_ca_ind <- readr::read_rds(
  here('data', 'cohort', 'ca_ind.rds')
)
dft_reg <- readr::read_rds(
  here('data', 'cohort', 'reg.rds')
)

# This file was created by going through the list of regimens at a meeting
#   one day.  The list was expanded on Feb 21 using the method of statistician
#   guesswork.
dft_neoadj <- readr::read_csv(
  here('data', 'drug_mapping', 'neoadjuvant_regimens_2024_02_21.csv'),
  show_col_types = F
)

dft_neoadj %<>% select(regimen_drugs, valid_neoadjuvant) # leaves off comments

dft_neoadj %<>% rename(valid_neoadjuvant_drug = valid_neoadjuvant) # for clarity

# Ensure all entries are 'yes' or 'no' only.
chk_neoadj_col <- dft_neoadj %>%
  pull(valid_neoadjuvant_drug) %>%
  `%in%`(., c("yes", "no")) %>%
  all(.)

if (!chk_neoadj_col) {
  cli::cli_abort(
    "Invalid entries in the input neoadjuvant CSV (valid_neoadjuvant_drug column)"
  )
}


# This neoadj drug table was created by hand.  Some of these regimens have since
#   been removed or corrected in the data.  The main example is regimens with
#.  both carbo and cis, but it turns out that one of those was a zero day use.
# We now correct the list to only include drugs in the data:
dft_neoadj %<>%
  filter(
    regimen_drugs %in% unique(pull(dft_reg, regimen_drugs))
  )

dft_neoadj %<>%
  arrange(regimen_drugs) %>%
  mutate(
    valid_neoadjuvant_drug = case_when(
      valid_neoadjuvant_drug %in% "yes" ~ T,
      valid_neoadjuvant_drug %in% "no" ~ F,
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
    valid_neoadjuvant_drug,
    dx_reg_start_int_yrs,
    dx_reg_end_all_int_yrs
  )

dft_dmet_timing <- get_dmet_time(dft_ca_ind) %>%
  select(record_id, ca_seq, dmet_time_yrs = dx_dmet_yrs)

dft_reg_neo <- left_join(
  dft_reg_neo,
  dft_dmet_timing,
  by = c("record_id", "ca_seq")
)

dft_reg_neo %<>%
  mutate(
    met_ever = if_else(is.na(dmet_time_yrs), F, T)
  )


dft_reg_neo %<>%
  filter(!str_detect(regimen_drugs, "Investigational Drug")) %>%
  mutate(
    tt_met_reg_start_yrs = dmet_time_yrs - dx_reg_start_int_yrs,
    tt_met_reg_end_yrs = dmet_time_yrs - dx_reg_end_all_int_yrs
  )

# We've updated our defintion.  Neoadjuvant is now defined by a drug use
#   which is on our list of regimens (platinum-based chemo) AND either:
#   1. drug use started more than 3 months before met diagnosis.
#   2. the participant never developed metastasis.
# There are some obvious problems with this, notably that the classification
#   of a drug as met or not could change if the person is diagnosed soon
#.  (in that three month no man's land).

dft_reg_neo %<>%
  mutate(
    is_neoadjuvant = case_when(
      !met_ever & valid_neoadjuvant_drug ~ T, # new case
      # small tolerance here, half a day, for rounding errors.
      valid_neoadjuvant_drug %in%
        T &
        tt_met_reg_start_yrs >= 0.25 - (0.5 / 365.25) ~
        T,
      T ~ F
    )
  )

dft_multi_neo <- dft_reg_neo %>%
  filter(is_neoadjuvant) %>%
  count(record_id, ca_seq, sort = T) %>%
  filter(n > 1)
# if (nrow(dft_multi_neo)) {
#   cli::cli_alert_info(
#     "Cases with multiple neoadjuvant regimens:"
#   )
#   print(dft_multi_neo)
# }

lev_neoadj <- c(
  "Neoadjuvant, started >=3 months prior to met",
  "Neoadjuvant, never-met patient",
  "Not neoadjuvant treatment"
)

dft_reg_neo %<>%
  mutate(
    neoadj_f = case_when(
      is.na(is_neoadjuvant) ~ lev_neoadj[3],
      !is_neoadjuvant ~ lev_neoadj[3],

      is_neoadjuvant & met_ever ~ lev_neoadj[1],
      is_neoadjuvant & !met_ever ~ lev_neoadj[2],
      T ~ NA_character_ # hopefully none.
    ),
    neoadj_f = factor(neoadj_f, levels = lev_neoadj)
  )
chk_met_neoadj_f <- dft_reg_neo %>%
  filter(is.na(neoadj_f)) %>%
  nrow(.)

if (chk_met_neoadj_f > 0) {
  cli::cli_abort(
    "There are regimens that were not categorized for 'met_neoadj_f' column."
  )
}


readr::write_rds(
  dft_reg_neo,
  file = here('data', 'dmet', 'neoadjuvant_status_reg.rds')
)


lev_neoadj_case <- lev_neoadj
lev_neoadj_case[3] <- "No neoadjuvant treatment"

dft_reg_neo_cases <- dft_reg_neo %>%
  group_by(record_id, ca_seq) %>%
  summarize(
    neoadj_case_f = case_when(
      any(neoadj_f %in% lev_neoadj[1]) ~ lev_neoadj_case[1],
      any(neoadj_f %in% lev_neoadj[2]) ~ lev_neoadj_case[2],
      T ~ lev_neoadj_case[3]
    ),
    tt_met_first_neoadj_start = suppressWarnings(
      min(tt_met_reg_start_yrs[is_neoadjuvant], na.rm = T)
    ),
    .groups = "drop"
  )

# Fix the +Inf returns from min above:
dft_reg_neo_cases %<>%
  mutate(
    tt_met_first_neoadj_start = if_else(
      is.infinite(tt_met_first_neoadj_start),
      NA_real_,
      tt_met_first_neoadj_start
    )
  )


# One more thing needs to be done here to make this more sensible:  We need to
#   add back in the people who had zero regimens.  Then they will be
#   correctly coded as having no neoadjuvant therapy.

dft_no_reg_addins <- dft_ca_ind %>%
  filter(!(record_id %in% dft_reg_neo_cases$record_id)) %>%
  select(record_id, ca_seq) %>%
  mutate(
    neoadj_case_f = lev_neoadj_case[3]
  )

dft_reg_neo_cases %<>%
  bind_rows(., dft_no_reg_addins) %>%
  mutate(
    neoadj_case_f = factor(neoadj_case_f, levels = lev_neoadj_case)
  )

dft_reg_neo_cases %<>%
  mutate(
    neoadj_case_binary = case_when(
      neoadj_case_f %in% lev_neoadj_case[1:2] ~ "Neoadj. exposure",
      T ~ "No Neoadj. exposure"
    )
  )

readr::write_rds(
  dft_reg_neo_cases,
  file = here('data', 'dmet', 'neoadjuvant_status_case.rds')
)
