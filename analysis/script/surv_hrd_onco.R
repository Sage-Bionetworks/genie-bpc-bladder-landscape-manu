# Create analysis dataset for a specific case identified by our physicians:
# Oncogenic HRD mutation or not is the main predictor (HRD = homologous repair deficiency)
# OS from first systemic therapy after metastasis is the main outcome.
# Drug groups: ABT + abi/enza and ABT + docetaxel  
# With and without risk set adjustment


library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source) # also load
dft_ca_ind <- readr::read_rds(here('data', 'cohort', "ca_ind.rds"))
dft_reg <- readr::read_rds(here('data', 'cohort', "reg.rds"))
dft_cpt <- readr::read_rds(here('data', 'cohort', "cpt_aug.rds"))
dft_alt <- readr::read_rds(here('data', 'genomic','alterations.rds'))



# Find the people who had a metastasis:
dft_met_timing <- get_dmet_time(dft_ca_ind)




dft_onco_hrd <- dft_alt %>% 
  filter(pathway_ddr %in% "HR") %>%
  filter(oncogenic %in% c("Likely Oncogenic", "Oncogenic"))

# Add in the relevant stuff from CPT data:
dft_onco_hrd <- dft_cpt %>%
  select(
    record_id, ca_seq,
    dx_path_proc_cpt_yrs, # interval from dx to pathology procedure of this CPT. 
    dx_cpt_rep_yrs, # interval from dx to report date of this CPT.
    cpt_genie_sample_id
  ) %>%
  left_join(
    dft_onco_hrd,
    .,
    by = c(sample_id = "cpt_genie_sample_id")
  )

# Find the people who had a metastasis and took a systemic therapy afterward:
dft_post_met_reg <- dft_reg %>%
  select(
    record_id, ca_seq, contains("regimen_number"),
    regimen_drugs,
    dx_reg_start_int_yrs,
    dx_reg_end_all_int,
    os_g_status, # g = regimen.  Why?  I just work here.
    tt_os_g_yrs
  ) %>%
  left_join(
    dft_met_timing,
    ., 
    by = c("record_id", "ca_seq")
  )


dft_post_met_reg %<>%
  # half day tolerance on the cutoff:
  filter(dx_reg_start_int_yrs >= (dx_dmet_yrs - 0.5 / 365.25))

dft_first_post_met_t <- dft_post_met_reg %>%
  group_by(record_id, ca_seq) %>%
  summarize(dx_first_post_met_reg_yrs = min(dx_reg_start_int_yrs, na.rm = T), .groups = "drop")

dft_first_cpt <- get_first_cpt(dft_ca_ind, dft_cpt) %>% rename(dx_first_cpt_rep_yrs = dx_cpt_rep_yrs)

dft_onco_hrd <- left_join(
  dft_onco_hrd,
  dft_first_post_met_t,
  by = c('record_id', 'ca_seq')
) %>%
  left_join(
    .,
    dft_first_cpt,
    by = c("record_id", 'ca_seq')
  )

dft_onco_hrd %<>% 
  mutate(
    # risk set entry happens when BOTH they have a CPT test and take a post-met drug.
    dx_entry = case_when(
      # they have to HAVE both those events, otherwise they did not enter:
      is.na(dx_first_cpt_rep_yrs) | is.na(dx_first_post_met_reg_yrs) ~ NA_real_,
      T ~ pmax(dx_first_cpt_rep_yrs, dx_first_post_met_reg_yrs, na.rm = T)
    )
  )

# dft_onco_hrd %>% 
#   ggplot(., aes(x = dx_path_proc_cpt_yrs, y = dx_entry)) + 
#   geom_point() + coord_equal()


readr::write_rds(
  dft_onco_hrd,
  here('data', 'outcome', 'surv_met_hrd', 'alt_onco_hrd.rds')
)





dft_onco_hrd_flags <- dft_onco_hrd %>%
  group_by(record_id, ca_seq) %>%
  summarize(
    hrd_before_pm_reg = sum(dx_path_proc_cpt_yrs < dx_first_post_met_reg_yrs, na.rm = T) >= 1,
    hrd_before_entry = sum(dx_path_proc_cpt_yrs < dx_entry, na.rm = T) >= 1,
    .groups = "drop"
  )

readr::write_rds(
  dft_onco_hrd_flags,
  here('data', 'outcome', 'surv_met_hrd', 'alt_onco_hrd_flags.rds')
)




dft_met_hrd_surv <- dft_post_met_reg %>%
  group_by(record_id, ca_seq) %>%
  arrange(dx_reg_start_int_yrs) %>%
  slice(1) %>%
  ungroup(.)

# temp_missing_cohort <- anti_join(
#   dft_onco_hrd_flags,
#   dft_met_hrd_surv,
#   by = c('record_id', 'ca_seq')
# )


dft_reg %>%
  filter(record_id %in% 'GENIE-DFCI-000893') %>% 
  glimpse

dft_met_hrd_surv <- left_join(
  dft_met_hrd_surv,
  dft_onco_hrd_flags,
  by = c("record_id", "ca_seq")
) 

dft_met_hrd_surv %<>%
  replace_na(
    list(
      hrd_before_pm_reg = 0,
      hrd_before_entry = 0
    )
  ) 

# Still need the cohort entry time for survival to work out:
dft_met_hrd_surv %<>%
  left_join(
    ., 
    dft_first_cpt,
    by = c('record_id', 'ca_seq')
  ) %>%
  # I just cant stand these names...
  rename(
    os_first_met_reg_status = os_g_status,
    tt_os_first_met_reg_yrs = tt_os_g_yrs
  ) %>%
  mutate(
    fmr_fcpt_yrs = dx_first_cpt_rep_yrs - dx_reg_start_int_yrs
  )






readr::write_rds(
  dft_met_hrd_surv,
  here('data', 'outcome', 'surv_met_hrd', 'met_hrd_surv.rds')
)




