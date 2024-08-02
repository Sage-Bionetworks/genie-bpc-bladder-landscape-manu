library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source) # also load

dft_met_ddr_surv <- readr::read_rds(
  here('data', 'survival', 'ddr_onco', 'met_ddr_surv.rds')
)

dir_out <- here('data', 'survival', 'ddr_onco')




surv_obj_os_fmr <- with(
  dft_met_ddr_surv,
  Surv(
    time = fmr_fcpt_yrs,
    time2 = tt_os_first_met_reg_yrs,
    event = os_first_met_reg_status
  )
)

# Some survival code I should probaly move later on:
dft_cox_univariate <- coxph(
  Surv(
    time = fmr_fcpt_yrs,
    time2 = tt_os_first_met_reg_yrs,
    event = os_first_met_reg_status
  ) ~ ddr_before_entry,
  data = dft_met_ddr_surv
) %>%
  broom::tidy(., conf.int = T)

readr::write_rds(
  dft_cox_univariate,
  here(dir_out, "cox_tidy_univariate.rds")
)








dft_met_ddr_surv_mod_ready <- dft_met_ddr_surv %>%
  mutate(dob_reg_start_yrs = dob_reg_start_days / 365.25) %>%
  select(
    fmr_fcpt_yrs,
    tt_os_first_met_reg_yrs,
    os_first_met_reg_status,
    ddr_before_entry,
    de_novo_met,
    dob_reg_start_yrs,
    md_ecog_imp_num,
    institution,
    race_eth,
    female
  ) %>%
  rename(age_reg_start = dob_reg_start_yrs) %>%
  fastDummies::dummy_cols(
    select_columns = c("institution", "race_eth"),
    remove_selected_columns = T, 
    remove_most_frequent_dummy = T
  ) # easier for me.

readr::write_rds(
  dft_met_ddr_surv_mod_ready,
  here(dir_out, "met_ddr_surv_mod_ready.rds")
)

# graph to share
md.pattern(dft_met_ddr_surv_mod_ready, rotate.names = T)

# Some problems with this as we have ECOG data that's about half missing.
dft_cox_complete_cases <- coxph(
  Surv(
    time = fmr_fcpt_yrs,
    time2 = tt_os_first_met_reg_yrs,
    event = os_first_met_reg_status
  ) ~ .,
  data = dft_met_ddr_surv_mod_ready
) %>%
  broom::tidy(., conf.int = T) 

readr::write_rds(
  dft_cox_complete_cases,
  here(dir_out, "cox_tidy_comp_cases.rds")
)







# Do the multiple imputation, using every other covariate to predict the
#   missing ecog scores.

blocks <- construct.blocks(
  formulas = list(md_ecog_imp_num ~ .)
)
pm <- make.predictorMatrix(
  data = dplyr::select(
    dft_met_ddr_surv_mod_ready, 
    -c(
      # just taking away the survival variables for predicting md_ecog_imp_sum
      fmr_fcpt_yrs,
      tt_os_first_met_reg_yrs,
      os_first_met_reg_status
    )),
    blocks = blocks
  )

mids_surv <- mice(
  data = dft_met_ddr_surv_mod_ready, 
  maxit = 5, m = 15, seed = 2341, 
  predictorMatrix = pm, 
  blocks = blocks,
  print = FALSE
)

gg_imp <- plot_gg_strip(mids_surv, var = "md_ecog_imp_num",
                        pt_size = 1)

readr::write_rds(
  x = gg_imp,
  file = here(dir_out, 'cox_imputation_plot_md_ecog.rds')
)

# Put the multiply imputed datasets into a list column
dat_res <- tibble(
  dat = mice::complete(mids_surv, "all")
)


# A bit of a pain here
mira_surv <- with(
  mids_surv,
  coxph(
    Surv(
      time = fmr_fcpt_yrs,
      time2 = tt_os_first_met_reg_yrs,
      event = os_first_met_reg_status
    ) ~ .,
    data = dft_met_ddr_surv_mod_ready
  )
)

dft_cox_mult_imp <- pool(mira_surv) %>% summary(conf.int = T) %>% select(term, estimate, conf.low = `2.5 %`, conf.high = `97.5 %`, p.value) %>%
  mutate(level = 0.95) #alpha level, since I blurred that out above.

readr::write_rds(
  dft_cox_mult_imp,
  here(dir_out, "cox_tidy_mult_imp.rds")
)

dft_cox_all_mods <- bind_rows(
  mutate(dft_cox_univariate, model = "uni"),
  mutate(dft_cox_complete_cases, model = "cc"),
  mutate(dft_cox_mult_imp, model = "mi")
)

  
  









# Obviously we want to replace this with imputation:
dft_met_ddr_surv_sub <- dft_met_ddr_surv_sub[complete.cases(dft_met_ddr_surv_sub),] 


cli_abort("Need to fix the UHN & race unknown issue")

obj_surv <- with(
  dft_met_ddr_surv_sub,
  Surv(
    time = fmr_fcpt_yrs,
    time2 = tt_os_first_met_reg_yrs,
    event = os_first_met_reg_status
  )
)

fit <- cv.glmnet(
  x = as.matrix(
    select(dft_met_ddr_surv_sub,
           -c(fmr_fcpt_yrs, tt_os_first_met_reg_yrs, os_first_met_reg_status)
    )
  ),
  y = obj_surv,
  family = "cox",
  nfolds = 5,
  type.measure = "C",
  alpha = 0.5
)

coef(fit, s = "lambda.1se")

# I think if we did a regularized fit with multiple imputation we'd be in better shape.

