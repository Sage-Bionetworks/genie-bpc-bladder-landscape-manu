# Yet another iteration of this question, rephrased by Michal and Neil for the
# 2024 ASCO urinary abstract deadline.

library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source) # also load

dft_pt <- readr::read_rds(here('data', 'cohort', "pt.rds"))
dft_ca_ind <- readr::read_rds(here('data', 'cohort', "ca_ind.rds"))
dft_reg <- readr::read_rds(here('data', 'cohort', "reg.rds"))
dft_cpt <- readr::read_rds(here('data', 'cohort', "cpt_aug.rds"))
dft_alt <- readr::read_rds(here('data', 'genomic','alterations.rds'))
dft_med_onc <- readr::read_rds(here('data', 'cohort', "med_onc.rds"))
dft_med_onc %<>%
  augment_med_onc_imputed_ecog(., add_numeric = T)

dir_out <- here('data', 'survival', 'asco_2024_analysis')



# Find the people who had a metastasis:
dft_met_timing <- get_dmet_time(dft_ca_ind)

# Provided by Michal Sternschuss from a trial:
custom_ddr_list <- c(
  'ERCC2', 'ERCC5', 
  'BRCA1', 'BRCA2', 'RECQL4', 'RAD51C', 'ATM', 
  'ATR', 'FANCC'
)

dft_onco_ddr <- dft_alt %>% 
  # This way gets everything in the Pearl list:
  # filter(!(pathway_ddr %in% "None")) %>%
  # Or we can always go with the custom pull:
  filter(hugo %in% custom_ddr_list) %>%
  filter(oncogenic %in% c("Likely Oncogenic", "Oncogenic"))

# First question:  How many people in the cohort have an oncogenic DDR mutation?
dft_cohort_ddr <- dft_onco_ddr %>% 
  mutate(feature = paste0(hugo, "_", str_sub(tolower(alt_type), 1, 3))) %>%
  select(sample_id, feature) %>%
  mutate(alt = 1)

dft_cohort_ddr <- dft_cpt %>%
  select(sample_id = cpt_genie_sample_id, 
         record_id,
         ca_seq) %>%
  left_join(., dft_cohort_ddr, by = c('sample_id')) %>%
  select(record_id, feature, alt) %>%
  filter(!is.na(feature))

skel_cohort_ddr <- tidyr::expand_grid(
  record_id = unique(dft_ca_ind$record_id),
  feature = unique(dft_cohort_ddr$feature)
)

dft_cohort_ddr <- left_join(
  skel_cohort_ddr,
  dft_cohort_ddr,
  by = c('record_id', 'feature')
) %>%
  replace_na(list(alt = 0))

dft_cohort_ddr <- dft_cohort_ddr %>%
  # summarize so that if any sample is positive for a feature, we count them.
  group_by(record_id, feature) %>%
  summarize(alt = any(alt %in% 1), .groups = 'drop') %>%
  pivot_wider(names_from = feature, values_from = alt)

# Add a count for any DDR mut.
dft_cohort_ddr %<>%
  mutate(
    any_ddr = rowSums(
      (dft_cohort_ddr %>% select(-record_id))
    ) > 0
  ) %>%
  relocate(any_ddr, .after = record_id)

readr::write_rds(
  dft_cohort_ddr,
  file = here(dir_out, "cohort_ddr.rds")
)





# Back to the main thread for the rest of the questions.
  
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

# Find the people who had a metastasis and took a platinum therapy after:
dft_post_met_plat <- dft_reg %>%
  # A change here, any platinum rather than just GemCis or GemCarbo.
  filter(
    str_detect(regimen_drugs, "Cisplatin|Carboplatin")
  ) %>% 
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

dft_post_met_plat %<>%
  # half day tolerance on the cutoff:
  filter(dx_reg_start_int_yrs >= (dx_dmet_yrs - 0.5 / 365.25))

dft_first_post_met_t <- dft_post_met_plat %>%
  group_by(record_id, ca_seq) %>%
  summarize(dx_first_post_met_reg_yrs = min(dx_reg_start_int_yrs, na.rm = T), .groups = "drop")

dft_first_cpt <- get_first_cpt(dft_ca_ind, dft_cpt) %>% 
  rename(dx_first_cpt_rep_yrs = dx_cpt_rep_yrs)

dft_onco_ddr <- left_join(
  dft_onco_ddr,
  dft_first_post_met_t,
  by = c('record_id', 'ca_seq')
) %>%
  left_join(
    .,
    dft_first_cpt,
    by = c("record_id", 'ca_seq')
  )

dft_onco_ddr %<>% 
  mutate(
    # risk set entry happens when BOTH they have a CPT test and take a post-met drug.
    dx_entry = case_when(
      # they have to HAVE both those events, otherwise they did not enter:
      is.na(dx_first_cpt_rep_yrs) | is.na(dx_first_post_met_reg_yrs) ~ NA_real_,
      T ~ pmax(dx_first_cpt_rep_yrs, dx_first_post_met_reg_yrs, na.rm = T)
    )
  )


dft_onco_ddr_flags <- dft_onco_ddr %>%
  group_by(record_id, ca_seq) %>%
  # half day tolerance for rounding aroudn each of these:
  summarize(
    ddr_before_pm_reg = sum(dx_path_proc_cpt_yrs < dx_first_post_met_reg_yrs + 0.5/365, na.rm = T) >= 1,
    ddr_before_entry = sum(dx_path_proc_cpt_yrs < dx_entry + 0.5/365 , na.rm = T) >= 1,
    .groups = "drop"
  )



# Option 2: To select only first line therapy:
dft_lot <- readr::read_rds(here('data', 'dmet', 'lines_of_therapy', 'lot.rds'))
dft_met_ddr_surv <- dft_lot %>%
  filter(
    line_therapy %in% 1,
    str_detect(regimen_drugs, "Cisplatin|Carboplatin")
  ) %>%
  select(record_id, ca_seq, regimen_number) %>%
  left_join(
    .,
    dft_post_met_plat,
    by = c('record_id', 'ca_seq', 'regimen_number')
  )



 
dft_met_ddr_surv <- left_join(
  dft_met_ddr_surv,
  dft_onco_ddr_flags,
  by = c("record_id", "ca_seq")
) 

dft_met_ddr_surv %<>%
  replace_na(
    list(
      ddr_before_pm_reg = 0,
      ddr_before_entry = 0
    )
  ) 

# Still need the cohort entry time for survival to work out:
dft_met_ddr_surv %<>%
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





# Adding a few things here which might be valuable as confounders.
# Even though this isn't currently used I'm leaving it in because it's where
#   we should be going.
dft_extra_var <- readr::read_rds(
  here('data', 'cohort', 'time_invariant_model_factors.rds')
) %>%
  filter(!is.na(dx_dmet_yrs)) %>%
  select(record_id, ca_seq, dx_dmet_yrs, dob_ca_dx_yrs, institution,
         birth_year, race_eth, female)

rc_vec <- c('record_id', 'ca_seq')
dft_reg_start_dob <- left_join(
  select(dft_met_ddr_surv, record_id, ca_seq, dx_reg_start_int_yrs),
  select(dft_ca_ind, record_id, ca_seq, dob_ca_dx_yrs),
  by = rc_vec
) %>%
  mutate(
    dob_reg_start_yrs = dx_reg_start_int_yrs + dob_ca_dx_yrs,
    dob_reg_start_days = dob_reg_start_yrs * 365.25
  ) %>%
  select(record_id, dob_reg_start_days) 

dft_extra_var %<>%
  left_join(., dft_reg_start_dob, by = 'record_id')

dft_latest_med_onc <- get_med_onc_by_timing(
  dat_med_onc = dft_med_onc,
  dat_cutoff = dft_reg_start_dob,
  var_cutoff = "dob_reg_start_days",
  remove_missing = T
)

dft_extra_var %<>%
  left_join(
    .,
    select(dft_latest_med_onc, record_id, md_ecog_imp_num),
    by = "record_id"
  )

dft_met_ddr_surv <- left_join(
  dft_met_ddr_surv,
  select(dft_extra_var, -dx_dmet_yrs), # no need to double up.
  by = rc_vec
)

# At this point I noticed one case where both carboplatin and cisplatin were documented.
# However, the carboplatin started and ended on the same day, so I feel confident
#   removing it.
# This is absolutely a hard code, and not durable to changes.
dft_met_ddr_surv %<>%
  mutate(
    regimen_drugs = case_when(
      record_id %in% 'GENIE-DFCI-006944' & ca_seq %in% 3 & regimen_number %in% 2 ~ "Cisplatin, Paclitaxel",
      T ~ regimen_drugs
    )
  )


dft_met_ddr_surv %<>% 
  remove_trunc_gte_event(
    trunc_var = 'fmr_fcpt_yrs',
    event_var = 'tt_os_first_met_reg_yrs'
  )

dft_met_ddr_surv %<>% mutate(fmr_fcpt_yrs = ifelse(fmr_fcpt_yrs < 0, 0, fmr_fcpt_yrs))

dft_met_ddr_surv_grouped <- bind_rows(
  (dft_met_ddr_surv %>% mutate(analysis_group = "any_plat")),
  (dft_met_ddr_surv %>% 
     filter(str_detect(regimen_drugs, "Cisplatin")) %>%
     mutate(analysis_group = "cisplatin_based")),
  (dft_met_ddr_surv %>% 
     filter(str_detect(regimen_drugs, "Carboplatin")) %>%
     mutate(analysis_group = "carboplatin_based"))
)

dft_met_ddr_surv_grouped %<>%
  nest(.by = analysis_group)



dft_met_ddr_surv_grouped %<>%
  mutate(
    cox_ddr = purrr::map(
      .x = data,
      .f = \(z) {
        coxph(
          data = z,
          with(z, Surv(
            time = fmr_fcpt_yrs,
            time2 = tt_os_first_met_reg_yrs,
            event = os_first_met_reg_status
          )) ~ ddr_before_entry
        ) %>%
          broom::tidy(., conf.int = T, exponentiate = F) # will exp later.
      }
    )
  )


dft_met_ddr_surv_grouped %<>%
  mutate(
    median_surv = purrr::map(
      .x = data,
      .f = \(z) {
        survfit(
          data = z,
          with(z, Surv(
            time = fmr_fcpt_yrs,
            time2 = tt_os_first_met_reg_yrs,
            event = os_first_met_reg_status
          )) ~ ddr_before_entry
        ) %>%
          summary %>%
          `$`(.,'table') %>%
          as_tibble(., rownames = 'group') %>%
          rename(lower = `0.95LCL`, upper = `0.95UCL`)
      }
    )
  )
  





# Print out these results (move later on):
dft_cohort_ddr %>% 
  select(-record_id) %>%
  gtsummary::tbl_summary(.)

dft_met_ddr_surv_grouped %>%
  select(analysis_group, cox_ddr) %>%
  unnest(cox_ddr) %>%
  # put the estimates on the HR scale (not log HR):
  mutate(
    across(
      .cols = c(estimate, conf.low, conf.high),
      .fns = exp
    )
  ) %>%
  select(1,3,6:8)

dft_met_ddr_surv_grouped %>%
  select(analysis_group, median_surv) %>%
  unnest(median_surv) %>%
  # put the estimates in months:
  mutate(
    across(
      .cols = c(median, lower, upper),
      .fns = \(z) z * 12.0148 # previously 12, makes almost no difference.
    )
  ) %>%
  # for now we just care about the medians:
  select(1,2, median, lower, upper)
  
dft_cohort_ddr %>% 
  select(-record_id) %>%
  gtsummary::tbl_summary(.)


dft_met_ddr_surv_grouped %>% 
  slice(2) %>%
  pull(data) %>%
  `[[`(.,1) %>%
  count(ddr_before_entry)


# It just never ends.
dft_ca_ind %>%
  filter(record_id %in% dft_met_timing$record_id) %>%
  select(age_dx, ca_type) %>%
  tbl_summary(.)

dft_pt %>%
  filter(record_id %in% dft_met_timing$record_id) %>%
  select(naaccr_sex_code) %>%
  tbl_summary(.)
