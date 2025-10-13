library(purrr)
library(fs)
library(here)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

read_wrap <- function(p) {
  read_rds(
    file = here("data", 'cohort', paste0(p, ".rds"))
  )
}

dft_pt <- read_wrap("pt")
dft_ca_ind <- read_wrap("ca_ind")
dft_regimens <- read_wrap("reg")

dft_drug <- dft_ca_ind %>%
  select(record_id, ca_seq) %>%
  inner_join(
    dft_regimens,
    .,
    by = c("record_id", "ca_seq")
  ) %>%
  select(
    record_id:redcap_ca_index,
    drugs_drug_1:dx_drug_end_or_lastadm_int_5
  ) %>%
  mutate(across(
    c(contains("drugs_drug_")),
    .fn = ~ vec_cast(.x, to = "character")
  )) %>%
  mutate(
    across(
      .cols = drugs_startdt_int_1:dx_drug_end_or_lastadm_int_5,
      .fns = as_double
    )
  )


# split off the character values, we'll do those separately
#   to avoid casting issues.
dft_drug_char <- dft_drug %>%
  select(record_id, regimen_number, ca_seq, drugs_drug_1:drugs_drug_5)

dft_drug %<>%
  select(-contains("drugs_drug_")) %>%
  pivot_longer(
    cols = drugs_startdt_int_1:dx_drug_end_or_lastadm_int_5,
    names_to = "var",
    values_to = "value"
  ) %>%
  # Because "drug_num" is confusing when we have "drugs_num" in
  #   the raw data.
  # This is only an id for drug within regimen within person.
  mutate(
    drug_id = readr::parse_number(var),
    var = stringr::str_replace(var, "_[0-5]", "")
  ) %>%
  pivot_wider(names_from = "var", values_from = "value") %>%
  select(record_id, contains("regimen_number"), drug_id, everything())


# Now do it again for the cluster of character columns.
#   All the previous ones were integer/double.
dft_drug_char %<>%
  pivot_longer(
    cols = drugs_drug_1:drugs_drug_5,
    names_to = "var",
    values_to = "value"
  ) %>%
  # Because "drug_num" is confusing when we have "drugs_num" in
  #   the raw data.
  # This is only an id for drug within regimen within person.
  mutate(
    drug_id = readr::parse_number(var),
    var = stringr::str_replace(var, "_[0-5]", "")
  ) %>%
  # because there is only one column we don't need to pivot here,
  #  just rename.
  select(record_id, regimen_number, drug_id, drug = value)

dft_drug <-
  left_join(
    dft_drug,
    dft_drug_char,
    by = c("record_id", "regimen_number", "drug_id")
  ) %>%
  relocate(drug, .before = drugs_startdt_int)

# empty rows here have no meaning - it's just regimens with less
#   than 5 drugs which is not suprising or interesting.
dft_drug %<>%
  filter(!is.na(drug))

readr::write_rds(
  x = dft_drug,
  file = here('data', 'cohort', "drug.rds")
)
