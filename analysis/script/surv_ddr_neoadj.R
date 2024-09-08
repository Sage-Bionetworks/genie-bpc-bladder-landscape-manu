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







dft_neoadj_ddr <- dft_ca_ind_ddr %>% 
  filter(ca_tx_pre_path_stage %in% "Yes", !(stage_dx_iv %in% "Stage IV")) 

dft_neoadj_ddr <- dft_neoadj_ddr %>%
  select(
    record_id,
    naaccr_path_t_cd,
    ca_path_t_stage,
    naaccr_clin_t_cd,
    ca_clin_t_stage,
    
    naaccr_path_n_cd,
    ca_path_n_stage,
    naaccr_clin_n_cd,
    ca_clin_n_stage
  ) %>%
  # Just make the unknown entries into actual NA's for ease:
  mutate(
    across(
      .cols = -record_id,
      .fns = \(z) if_else(
        z %in% "Unknown",
        NA_character_,
        z
      )
    )
  )

dft_neoadj_ddr_long_before <- dft_neoadj_ddr %>%
  pivot_longer(
    cols = -record_id,
    values_to = "before"
  )
# Checking that worked:
purrr::walk(
  .x = names(dft_neoadj_ddr)[-1],
  .f = \(z) count(dft_neoadj_ddr, .data[[z]], sort = T) %>% print
)

# Remove all the redundant (with col names), inconsistent and troublingly error prone stuff
#   at the beginning.  We just want those digits.
dft_neoadj_ddr %<>%
  mutate(
    across(
      .cols = -record_id,
      .fns = \(z) 
      str_replace_all(
        z, "^[TtPpNnCc]*", ""
      )
    )
  )

# A few special cases:
dft_neoadj_ddr %<>%
  mutate(
    across(
      .cols = -record_id,
      .fns = \(z) {
        case_when(
          z %in% "IS" ~ "0",
          z %in% "X" ~ NA_character_,
          T ~ z
        )
      }
    )
  )

# Now let's get the number.  I'm expecting this to be the first char.
dft_neoadj_ddr %<>%
  mutate(
    across(
      .cols = -record_id,
      .fns = \(z) {
        str_extract(z, "^[0-9]")
      }
    )
  )

dft_neoadj_ddr_long_after <- dft_neoadj_ddr %>%
  pivot_longer(
    cols = -record_id,
    values_to = "after"
  )

# Check the mapping if desired:
# full_join(
#   dft_neoadj_ddr_long_before,
#   dft_neoadj_ddr_long_after,
#   by = c('record_id', 'name')
# ) %>%
#   select(name, before, after) %>%
#   distinct(.) %>%
#   print(n = 500)

dft_neoadj_ddr %<>%
  mutate(
    across(
      .cols = -record_id,
      .fns = as.numeric
    )
  )
  
dft_neoadj_ddr %<>%
  mutate(
    path_comb_n = case_when(
      !is.na(naaccr_path_n_cd) ~ naaccr_path_n_cd,
      T ~ ca_path_n_stage
    ),
    clin_comb_n = case_when(
      !is.na(naaccr_clin_n_cd) ~ naaccr_clin_n_cd,
      T ~ ca_clin_n_stage
    ),
    path_comb_t = case_when(
      !is.na(naaccr_path_t_cd) ~ naaccr_path_t_cd,
      T ~ ca_path_t_stage
    ),
    clin_comb_t = case_when(
      !is.na(naaccr_clin_t_cd) ~ naaccr_clin_t_cd,
      T ~ ca_clin_t_stage
    )
  )

dft_neoadj_ddr %<>%
  select(
    record_id,
    contains("comb")
  ) %>%
  mutate(
    diff_n = path_comb_n - clin_comb_n,
    diff_t = path_comb_t - clin_comb_t
  )



readr::write_rds(
  dft_neoadj_ddr,
  here(dir_output, 'ddr_neoadj_stage_diffs.rds')
)