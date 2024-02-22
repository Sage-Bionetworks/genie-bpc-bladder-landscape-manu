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
  "Mitomycin"
)
# Report those out so they can go in a report section:
readr::write_rds(
  x = excluded_drug_strings,
  file = here(
    'data', 'dmet', 'lines_of_therapy',
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

# Temporary note to help me find PD1 and PDL1 inhibitors.
dft_reg_dmet_s %>% 
  filter(str_detect(tolower(regimen_drugs), paste(
    c("pembrolizumab",
      "nivolumab",
      "cemiplimab",
      "dostarlimab",
      "retifanlimab",
      "toripalimab",
      "atezolizumab",
      "avelumab",
      "durvalumab"
    ),
    collapse = "|"
  ))) %>%
  tabyl(regimen_drugs)


# Variation classes are (in this context) drug regimens that would not constitute
#   a new regimen.  For example, carboplatin + gemcitabine is in the same class
#   as cisplatin + gemcitabine.  So if we see GemCis after GemCarbo then 
#   we don't count GemCis as a second line of therapy, it's a variation on the first.
dft_variation_classes <- tribble(
  ~var_class, ~regimen_drugs,
  # Class 1:  Platinum chemo with gemcitabine
  1, 'Cisplatin, Gemcitabine Hydrochloride',
  1, 'Carboplatin, Gemcitabine Hydrochloride',
  # Class 2: PD1 and PDL1 inhibitors without chemo
  2, "Atezolizumab", 
  2, "Avelumab",
  2, "Durvalumab",
  2, "Enfortumab Vedotin, Pembrolizumab",
  2, "Ipilimumab, Nivolumab",
  2, "Nivolumab",
  2, "Pembrolizumab"
)

if (!all( dft_variation_classes$regimen_drugs %in% dft_reg_dmet_s$regimen_drugs) ) {
  cli_abort("Variation class declaration error - 1+ not found in regimen data (typo likely)")
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
  here('data', 'dmet' , 'lines_of_therapy', 'lot.rds')
)


