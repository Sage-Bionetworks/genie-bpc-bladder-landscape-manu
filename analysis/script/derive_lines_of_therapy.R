# Description:  Derive lines of therapy starting from the regimen data on/after
#   metastatic diagnosis.

library(fs)
library(purrr)
library(here)
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

# Update from Dec 2023:  Keep the investigational drugs but label them
#   ALL as "investigational regimen"
dft_reg_dmet_s %<>%
  mutate(
    regimen_drugs = if_else(
      str_detect(regimen_drugs, "Investigational"),
      "Investigational Regimen",
      regimen_drugs
    )
  )

excluded_drug_strings <- c(
  "BCG Vaccine",
  "BCG Solution",
  "Mitomycin",
  # discussed adding this on Jan 27 meeting:
  "Autologous Melanoma Lysate Pulsed Dendritic Cell Vaccine, Other NOS"
)


# Report those out so they can go in a report section:
readr::write_rds(
  x = excluded_drug_strings,
  file = here(
    'data',
    'dmet',
    'lines_of_therapy',
    'excluded_drugs.rds'
  )
)

# Some caution is needed here.  This will exclude any regimen that contains those
#   drug strings at all.  This might be inadvisable if the regimen has other
#   valid drug entries.  As it stands right now, excluded regimens are mostly
#   monotherapies, so this works fine.
dft_reg_dmet_s %<>%
  filter(
    !str_detect(
      regimen_drugs,
      paste0(excluded_drug_strings, collapse = "|")
    )
  )

# We will also exclude those with any non-standard administration method.
# Currently this is empirically limited to intravesicular administration,
#   which is most common with the BCG vaccine, so this is almost fully redundant
#   with the above.
dft_reg_dmet_s %<>%
  filter(
    is.na(drugs_admin) # NA = standard in PRISSMM.
  )

# Variation classes are (in this context) drug regimens that would not constitute
#   a new regimen.  For example, carboplatin + gemcitabine is in the same class
#   as cisplatin + gemcitabine.  So if we see GemCis after GemCarbo then
#   we don't count GemCis as a second line of therapy, it's a variation on the first.
dft_variation_classes <- tribble(
  ~var_class,
  ~regimen_drugs,
  # Class 1:  Platinum chemo with gemcitabine
  1,
  'Cisplatin, Gemcitabine Hydrochloride',
  1,
  'Carboplatin, Gemcitabine Hydrochloride',
  # Class 2: PD1 and PDL1 inhibitors without chemo
  2,
  "Atezolizumab",
  2,
  "Avelumab",
  2,
  "Durvalumab",
  2,
  "Enfortumab Vedotin, Pembrolizumab",
  2,
  "Ipilimumab, Nivolumab",
  2,
  "Nivolumab",
  2,
  "Pembrolizumab",
  3,
  "Carboplatin, Gemcitabine Hydrochloride, Paclitaxel",
  3,
  "Cisplatin, Gemcitabine Hydrochloride, Paclitaxel",
  4,
  "Cisplatin, Etoposide",
  4,
  "Carboplatin, Etoposide"
)

if (
  !all(dft_variation_classes$regimen_drugs %in% dft_reg_dmet_s$regimen_drugs)
) {
  cli_abort(
    "Variation class declaration error - 1+ not found in regimen data (typo likely)"
  )
}


readr::write_rds(
  x = dft_variation_classes,
  file = here('data', 'dmet', 'lines_of_therapy', 'variation_classes.rds')
)

dft_reg_dmet_s %<>%
  left_join(
    .,
    dft_variation_classes,
    by = 'regimen_drugs'
  )

dft_lot <- dft_reg_dmet_s %<>%
  assign_lot(.)

readr::write_rds(
  x = dft_lot,
  here('data', 'dmet', 'lines_of_therapy', 'lot.rds')
)


# Addon for the moment:  Immunotherapy lines
io_lot <- dft_lot %>%
  filter(!is.na(line_therapy)) %>%
  mutate(
    has_atezolizumab = str_detect(tolower(regimen_drugs), 'atezolizumab'),
    has_pembrolizumab = str_detect(tolower(regimen_drugs), 'pembrolizumab'),
    has_avelumab = str_detect(tolower(regimen_drugs), 'avelumab'),
    has_durvalumab = str_detect(tolower(regimen_drugs), 'durvalumab'),
    has_ipilizumab = str_detect(tolower(regimen_drugs), 'ipilizumab'),
    has_nivolumab = str_detect(tolower(regimen_drugs), 'nivolumab'),
    has_any_io = has_atezolizumab |
      has_pembrolizumab |
      has_avelumab |
      has_durvalumab |
      has_ipilizumab |
      has_nivolumab
  ) %>%
  mutate(
    line_therapy = if_else(line_therapy > 3, "4+", as.character(line_therapy))
  ) %>%
  group_by(line_therapy) %>%
  summarize(
    across(
      .cols = matches('has_'),
      .fns = sum
    )
  )

readr::write_csv(
  io_lot,
  here('analysis', 'explore', 'io_lot.csv')
)
