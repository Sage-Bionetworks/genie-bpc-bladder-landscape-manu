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

# Variation classes are (in this context) drug regimens that would not constitute
#   a new regimen.  For example, carboplatin + gemcitabine is in the same class
#   as cisplatin + gemcitabine.  So if we see GemCis after GemCarbo then 
#   we don't count GemCis as a second line of therapy, it's a variation on the first.
dft_variation_classes <- tribble(
  ~var_class, ~regimen_drugs,
  1, 'Cisplatin, Gemcitabine Hydrochloride',
  1, 'Carboplatin, Gemcitabine Hydrochloride'
)

dft_reg_dmet_s %<>%
  left_join(
    .,
    dft_variation_classes,
    by = 'regimen_drugs'
  )

# It is assumed this has already been given a column 'var_classes'.
assign_lot <- function(
    dat_reg
) {
  if (!('var_class' %in% colnames(dft_reg_dmet_s))) {
    dat_reg %<>% 
      mutate(var_class = NA_character_)
  }
  
  dat_reg %<>%
    arrange(record_id, ca_seq, regimen_number)
  
  dat_reg %<>%
    group_by(record_id, ca_seq) %>%
    mutate(
      line_therapy = lot_helper_one_case(
        var_class = var_class
      )
    ) %>%
    ungroup(.)
  
  return(dat_reg)
  
}

dft_reg_dmet_s %>%
  assign_lot(.) %>%
  filter(record_id %in% "GENIE-DFCI-003553") %>%
  select(regimen_drugs, line_therapy)




lot_helper_one_case <- function(var_class) {
  lot <- c()
  observed_var_classes <- c()
  lot_iter <- 0L
  # for loop horror:
  for (i in 1:length(var_class)) {
    if (var_class[i] %in% observed_var_classes) {
      lot[i] <- NA_real_
    } else {
      lot_iter <- lot_iter + 1L
      lot[i] <- lot_iter
      if (!is.na(var_class[i])) { 
        observed_var_classes <- c(observed_var_classes, var_class[i])
      }
    }
  }
  return(lot)
}

dft_reg_dmet_s %>% 
  filter(record_id %in% "GENIE-DFCI-003553") %>% 
  # mutate(var_class = NA_character_) %>%
  pull(var_class) %>%
  lot_helper_one_case(.)


# Todo: Declare equivalent classes of drugs, maybe start with GemCis and GemCarbo.