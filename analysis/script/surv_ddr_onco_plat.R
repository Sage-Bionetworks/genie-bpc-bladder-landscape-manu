# Create analysis dataset for a specific case identified by our physicians:
# Oncogenic DDR from first platinum therapy after metastasis is the main outcome.

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

dir_out <- here('data', 'survival', 'ddr_onco')



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
  filter(regimen_drugs %in% "Cisplatin, Gemcitabine Hydrochloride") %>%
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

# dft_onco_ddr %>% 
#   ggplot(., aes(x = dx_path_proc_cpt_yrs, y = dx_entry)) + 
#   geom_point() + coord_equal()


readr::write_rds(
  dft_onco_ddr,
  here(dir_out, 'alt_onco_ddr.rds')
)





dft_onco_ddr_flags <- dft_onco_ddr %>%
  group_by(record_id, ca_seq) %>%
  summarize(
    ddr_before_pm_reg = sum(dx_path_proc_cpt_yrs < dx_first_post_met_reg_yrs, na.rm = T) >= 1,
    ddr_before_entry = sum(dx_path_proc_cpt_yrs < dx_entry, na.rm = T) >= 1,
    .groups = "drop"
  )

readr::write_rds(
  dft_onco_ddr_flags,
  here(dir_out, 'alt_onco_ddr_flags.rds')
)


# Two Options here.

# Option 1: To just select the first regimen, no matter how it falls in lines:
# dft_met_ddr_surv <- dft_post_met_plat %>%
#   group_by(record_id, ca_seq) %>%
#   arrange(dx_reg_start_int_yrs) %>%
#   slice(1) %>%
#   ungroup(.)

# Option 2: To select only first line therapy:
dft_lot <- readr::read_rds(here('data', 'dmet', 'lines_of_therapy', 'lot.rds'))
dft_met_ddr_surv <- dft_lot %>%
  filter(
    line_therapy %in% 1,
    regimen_drugs %in% "Cisplatin, Gemcitabine Hydrochloride"
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





# Adding a few things here which might be valuable as confounders:
dft_extra_var <- readr::read_rds(
  here('data', 'cohort', 'time_invariant_model_factors.rds')
) %>%
  filter(is.na(dx_dmet_yrs)) %>%
  select(record_id, ca_seq, dx_dmet_yrs, dob_ca_dx_yrs, institution,
         birth_year, race_eth, female)

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


readr::write_rds(
  dft_met_ddr_surv,
  here(dir_out, 'met_ddr_surv.rds')
)








# Create the basic survival plots

dft_met_ddr_surv %<>% 
  remove_trunc_gte_event(
    trunc_var = 'fmr_fcpt_yrs',
    event_var = 'tt_os_first_met_reg_yrs'
  )

dft_met_ddr_surv %<>% mutate(fmr_fcpt_yrs = ifelse(fmr_fcpt_yrs < 0, 0, fmr_fcpt_yrs))

surv_obj_os_fmr <- with(
  dft_met_ddr_surv,
  Surv(
    time = fmr_fcpt_yrs,
    time2 = tt_os_first_met_reg_yrs,
    event = os_first_met_reg_status
  )
)

dft_met_ddr_surv %<>%
  mutate(
    ddr_disp = case_when(
      ddr_before_entry ~ "Oncogenic DDR",
      T ~ "No Onco. DDR"
    )
  )

gg_os_fmr_ddr <- plot_one_survfit(
  dat = dft_met_ddr_surv,
  surv_form = surv_obj_os_fmr ~ ddr_disp,
  plot_title = "OS from first line (metastatic) GemCis",
  plot_subtitle = "Adjusted for (independent) delayed entry"
)

readr::write_rds(
  gg_os_fmr_ddr,
  file = here(dir_out, "gg_met_ddr.rds")
)


# Make a truncated version:

gg_os_fmr_ddr_aacr_ss24 <- plot_one_survfit(
  dat = dft_met_ddr_surv,
  surv_form = surv_obj_os_fmr ~ ddr_disp,
  plot_title = "OS from first line (metastatic) GemCis",
  plot_subtitle = "Adjusted for (independent) delayed entry",
  x_breaks = seq(0, 100, by = 0.5)
) + 
  coord_cartesian(xlim = c(0,5))

ggsave(
  plot = gg_os_fmr_ddr_aacr_ss24,
  height = 4, width = 7,
  filename = here('output', 'aacr_ss24', 'img', '03_met_ddr.pdf')
)









# Had a question later on about how many people have gemcis as second line therapy.
# We can count this up:
dft_met_ddr_surv_2nd_line <- dft_lot %>%
  filter(
    line_therapy %in% 2, # only chage here - rest is copy-paste
    regimen_drugs %in% "Cisplatin, Gemcitabine Hydrochloride"
  ) %>%
  select(record_id, ca_seq, regimen_number) %>%
  left_join(
    .,
    dft_post_met_plat,
    by = c('record_id', 'ca_seq', 'regimen_number')
  ) %>%
  left_join(
    .,
    dft_onco_ddr_flags,
    by = c("record_id", "ca_seq")
  ) 
dft_met_ddr_surv_2nd_line %>% count(ddr_before_pm_reg) # NA -> false in later step.
# I'm just going to add this as a manual note rather than saving this dataset.




