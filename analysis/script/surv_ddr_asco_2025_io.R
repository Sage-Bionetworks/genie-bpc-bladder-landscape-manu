# Yet another iteration of this question, rephrased by Michal and Neil for the
# 2024 ASCO urinary abstract deadline.

library(purrr)
library(here)
library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source) # also load

dft_pt <- readr::read_rds(here('data', 'cohort', "pt.rds"))
dft_ca_ind <- readr::read_rds(here('data', 'cohort', "ca_ind.rds"))
dft_reg <- readr::read_rds(here('data', 'cohort', "reg.rds"))
dft_cpt <- readr::read_rds(here('data', 'cohort', "cpt_aug.rds"))
dft_alt <- readr::read_rds(here('data', 'genomic', 'alterations.rds'))
dft_med_onc <- readr::read_rds(here('data', 'cohort', "med_onc.rds"))
dft_med_onc %<>%
  augment_med_onc_imputed_ecog(., add_numeric = T)
dft_lot <- readr::read_rds(here('data', 'dmet', 'lines_of_therapy', 'lot.rds'))


dir_out <- here('data', 'survival', 'asco_2024_analysis')


custom_ddr_list <- c(
  'ERCC2',
  'ERCC5',
  'BRCA1',
  'BRCA2',
  'RECQL4',
  'RAD51C',
  'ATM',
  'ATR',
  'FANCC'
)
dft_onco_ddr <- dft_alt %>%
  filter(hugo %in% custom_ddr_list) %>%
  filter(oncogenic %in% c("Likely Oncogenic", "Oncogenic"))

dft_cohort_ddr <- readr::read_rds(
  file = here(dir_out, "cohort_ddr.rds")
)


# Back to the main thread for the rest of the questions.

# Add in the relevant stuff from CPT data:
dft_onco_ddr <- dft_cpt %>%
  select(
    record_id, #ca_seq,
    dob_cpt_report_days,
    path_proc_cpt_rep_days,
    # dx_path_proc_cpt_yrs, # interval from dx to pathology procedure of this CPT.
    # dx_cpt_rep_yrs, # interval from dx to report date of this CPT.
    cpt_genie_sample_id
  ) %>%
  left_join(
    dft_onco_ddr,
    .,
    by = c(sample_id = "cpt_genie_sample_id")
  )

dft_met_timing <- get_dmet_time(dft_ca_ind)


dft_relevant_reg <- dft_lot %>%
  filter(
    str_detect(regimen_drugs, "Atezo|Pembro"),
    line_therapy %in% 2
  ) %>%
  left_join(
    .,
    select(
      dft_reg,
      record_id,
      ca_seq,
      regimen_number,
      dx_reg_start_int_yrs,
      os_g_status,
      tt_os_g_yrs
    ),
    by = c("record_id", 'ca_seq', 'regimen_number')
  )

dft_first_cpt <- get_first_cpt(dft_ca_ind, dft_cpt) %>%
  rename(dx_first_cpt_rep_yrs = dx_cpt_rep_yrs)

first_ddr_alt <- dft_onco_ddr %>%
  group_by(record_id) %>%
  arrange(dob_cpt_report_days) %>%
  slice(1) %>%
  ungroup(.) %>%
  select(
    record_id,
    dob_ddr_alt = dob_cpt_report_days,
    path_proc_cpt_rep_days
  )


dft_relevant_reg <- left_join(
  dft_relevant_reg,
  select(dft_ca_ind, record_id, ca_seq, dob_ca_dx_days),
  by = c('record_id', 'ca_seq')
) %>%
  left_join(
    .,
    first_ddr_alt,
    by = 'record_id'
  )

dft_relevant_reg %<>%
  mutate(
    dx_ddr_alt = dob_ddr_alt - dob_ca_dx_days,
    dx_ddr_alt_yrs = dx_ddr_alt / 365.25,
    dx_path_proc_ddr_alt = dx_ddr_alt - path_proc_cpt_rep_days,
    dx_path_proc_ddr_alt_yrs = dx_path_proc_ddr_alt / 365.25
  ) %>%
  select(-c(dob_ca_dx_days, dob_ddr_alt, dx_ddr_alt, dx_path_proc_ddr_alt))


dft_relevant_reg %<>%
  left_join(
    .,
    dft_first_cpt,
    by = c("record_id", 'ca_seq')
  )


dft_relevant_reg %<>%
  mutate(
    # risk set entry happens when BOTH they have a CPT test and take a post-met drug.
    dx_entry = case_when(
      # they have to HAVE both those events, otherwise they did not enter:
      is.na(dx_first_cpt_rep_yrs) | is.na(dx_reg_start_int_yrs) ~ NA_real_,
      T ~ pmax(dx_first_cpt_rep_yrs, dx_reg_start_int_yrs, na.rm = T)
    )
  )

# ggplot(
#   data = dft_relevant_reg,
#   aes(x = dx_ddr_alt_yrs, y = dx_entry)
# ) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0)

dft_relevant_reg <- dft_relevant_reg %>%
  mutate(
    ddr_alt = !is.na(dx_ddr_alt_yrs),
    ddr_before_pm_reg = case_when(
      !ddr_alt ~ F,
      dx_ddr_alt_yrs < dx_reg_start_int_yrs + 0.5 / 365 ~ T,
      T ~ F
    ),
    ddr_before_entry = case_when(
      !ddr_alt ~ F,
      dx_ddr_alt_yrs < dx_entry + 0.5 / 365 ~ T,
      T ~ F
    )
  ) %>%
  select(-ddr_alt)


dft_met_ddr_surv <- dft_relevant_reg

dft_met_ddr_surv %<>%
  replace_na(
    list(
      ddr_before_pm_reg = 0,
      ddr_before_entry = 0
    )
  )


dft_met_ddr_surv %<>%
  rename(
    os_reg_status = os_g_status,
    tt_os_rel_reg_yrs = tt_os_g_yrs
  ) %>%
  mutate(
    reg_fcpt_yrs = dx_first_cpt_rep_yrs - dx_reg_start_int_yrs
  )


# Adding a few things here which might be valuable as confounders.
# Even though this isn't currently used I'm leaving it in because it's where
#   we should be going.
dft_extra_var <- readr::read_rds(
  here('data', 'cohort', 'time_invariant_model_factors.rds')
) %>%
  filter(!is.na(dx_dmet_yrs)) %>%
  select(
    record_id,
    ca_seq,
    dx_dmet_yrs,
    dob_ca_dx_yrs,
    institution,
    birth_year,
    race_eth,
    female
  )

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


dft_met_ddr_surv %<>%
  remove_trunc_gte_event(
    trunc_var = 'reg_fcpt_yrs',
    event_var = 'tt_os_rel_reg_yrs'
  )

dft_met_ddr_surv %<>%
  mutate(reg_fcpt_yrs = ifelse(reg_fcpt_yrs < 0, 0, reg_fcpt_yrs))

dft_met_ddr_surv_grouped <- bind_rows(
  (dft_met_ddr_surv %>% mutate(analysis_group = "any_io")),
  (dft_met_ddr_surv %>%
    filter(str_detect(regimen_drugs, "Pembrolizumab")) %>%
    mutate(analysis_group = "pembro_based")),
  (dft_met_ddr_surv %>%
    filter(str_detect(regimen_drugs, "Atezolizumab")) %>%
    mutate(analysis_group = "atezo_based"))
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
          with(
            z,
            Surv(
              time = reg_fcpt_yrs,
              time2 = tt_os_rel_reg_yrs,
              event = os_reg_status
            )
          ) ~
            ddr_before_entry
        ) %>%
          broom::tidy(., conf.int = T, exponentiate = F) # will exp later.
      }
    )
  )


# Will come back to this:
dft_met_ddr_surv_grouped %<>%
  mutate(
    median_surv = purrr::map(
      .x = data,
      .f = \(z) {
        survfit(
          data = z,
          with(
            z,
            Surv(
              time = reg_fcpt_yrs,
              time2 = tt_os_rel_reg_yrs,
              event = os_reg_status
            )
          ) ~
            ddr_before_entry
        ) %>%
          summary %>%
          `$`(., 'table') %>%
          as_tibble(., rownames = 'group') %>%
          rename(lower = `0.95LCL`, upper = `0.95UCL`)
      }
    )
  )


# # Print out these results (move later on):
dft_cohort_ddr %>%
  select(-record_id) %>%
  gtsummary::tbl_summary(.)

dft_met_ddr_surv_grouped %>%
  mutate(
    pt_counts = purrr::map(
      .x = data,
      .f = \(x) count(x, ddr_before_entry)
    )
  ) %>%
  select(analysis_group, pt_counts) %>%
  unnest(pt_counts)

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
  select(1, 3, 6:8)

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
  select(1, 2, median, lower, upper)

dft_cohort_ddr %>%
  select(-record_id) %>%
  gtsummary::tbl_summary(.)


# dft_met_ddr_surv_grouped %>%
#   slice(2) %>%
#   pull(data) %>%
#   `[[`(.,1) %>%
#   count(ddr_before_entry)

# It just never ends - a few more asks here:
dft_ca_ind %>%
  filter(record_id %in% dft_met_timing$record_id) %>%
  select(age_dx, ca_type) %>%
  tbl_summary(.)

dft_pt %>%
  filter(record_id %in% dft_met_timing$record_id) %>%
  select(naaccr_sex_code) %>%
  tbl_summary(.)


surv_obj_os_fmr <- with(
  dft_met_ddr_surv,
  Surv(
    time = reg_fcpt_yrs,
    time2 = tt_os_rel_reg_yrs,
    event = os_reg_status
  )
)

dft_met_ddr_surv %<>%
  mutate(
    ddr_disp = case_when(
      ddr_before_entry ~ "Oncogenic DDR",
      T ~ "No Onco. DDR"
    )
  )

gg_os_2l_io_ddr <- plot_one_survfit(
  dat = dft_met_ddr_surv,
  surv_form = surv_obj_os_fmr ~ ddr_disp,
  plot_title = "OS from second line pembrolizumab or atezolizumab",
  plot_subtitle = "Adjusted for (independent) delayed entry",
  x_breaks = seq(0, 10, by = 0.5),
  x_exp = 0.05
)

# A no-risktable version for the poster.
gg_os_2l_io_ddr_no_rt <- plot_one_survfit_no_risktable(
  dat = dft_met_ddr_surv,
  surv_form = surv_obj_os_fmr ~ ddr_disp,
  plot_title = "OS from second line pembrolizumab or atezolizumab",
  plot_subtitle = "Adjusted for (independent) delayed entry",
  x_breaks = seq(0, 10, by = 0.5),
  x_exp = 0.05
)


readr::write_rds(
  gg_os_2l_io_ddr,
  here('data', 'survival', 'ddr_onco', 'gg_os_2l_io_ddr.rds')
)

readr::write_rds(
  gg_os_2l_io_ddr_no_rt,
  here('data', 'survival', 'ddr_onco', 'gg_os_2l_io_ddr_no_rt.rds')
)
