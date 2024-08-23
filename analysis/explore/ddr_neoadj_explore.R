# In an August 2024 meeting we discussed whether we could add to the analysis
#   of DDR mutations by looking at the difference between clinical and path
#   staging in those individuals.


library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

dir_output <- here('data', 'survival', 'ddr_neoadj')

dft_ca_ind <- readr::read_rds(here('data', 'cohort', "ca_ind.rds"))
dft_reg <- readr::read_rds(here('data', 'cohort', "reg.rds"))
dft_lot <- readr::read_rds(here('data', 'dmet', 'lines_of_therapy', 'lot.rds'))
dft_cpt <- readr::read_rds(here('data', 'cohort', "cpt_aug.rds"))





# First step:  how many people we got with a ddr mutation?

custom_ddr_list <- c(
  'ERCC2', 'ERCC5', 
  'BRCA1', 'BRCA2', 'RECQL4', 'RAD51C', 'ATM', 
  'ATR', 'FANCC'
)

dft_onco_ddr <- dft_alt %>% 
  filter(hugo %in% custom_ddr_list) %>%
  filter(oncogenic %in% c("Likely Oncogenic", "Oncogenic"))

# Add in the relevant stuff from CPT data:
dft_onco_ddr <- dft_cpt %>%
  select(
    record_id, ca_seq,
    dx_path_proc_cpt_yrs, # interval from dx to pathology procedure of this CPT. 
    dx_cpt_rep_yrs, # interval from dx to report date of this CPT.
    cpt_genie_sample_id
  ) %>%
  left_join(
    dft_onco_ddr,
    .,
    by = c(sample_id = "cpt_genie_sample_id")
  )

# First step part 2: How many of those people have neoadjuvant therapy according
#   to the index cancer dataset?

dft_ca_ind_ddr <- dft_onco_ddr %>%
  select(record_id, ca_seq) %>%
  distinct(.) %>%
  left_join(., dft_ca_ind, by = c('record_id', 'ca_seq'))

readr::write_rds(
  dft_ca_ind_ddr,
  here(dir_output, 'ca_ind_ddr.rds')
)


