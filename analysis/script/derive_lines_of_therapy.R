# Description:  Derive lines of therapy starting from the regimen data on/after
#   metastatic diagnosis.

library(fs); library(purrr); library(here);
purrr::walk(.x = fs::dir_ls('R'), .f = source)

dft_reg_dmet_s <- readr::read_rds(
  here('data', 'dmet', 'reg_start_gte_dmet.rds')
)
dft_neoadj_case <- readr::read_rds(
  here('data', 'dmet', 'neoadjuvant_status_case.rds')
)

dft_reg_dmet_s %<>% 
  left_join(
    .,
    dft_neoadj_case,
    by = c("record_id", "ca_seq")
  )

dft_reg_dmet_s %<>%
  filter(!str_detect(regimen_drugs, "Investigational Drug"))

excluded_drug_strings <- c(
  "BCG Vaccine",
  "BCG Solution"
)
# Report those out so they can go in a report section:
readr::write_rds(
  x = excluded_drug_strings,
  file = here(
    'data', 'dmet', 'lines_of_therapy',
    'excluded_drugs.rds'
  )
)

dft_reg_dmet_s %<>%
  filter(
    !str_detect(
      regimen_drugs,
      paste0(excluded_drug_strings, collapse = "|")
    )
  )


# Todo: Declare equivalent classes of drugs, maybe start with GemCis and GemCarbo.