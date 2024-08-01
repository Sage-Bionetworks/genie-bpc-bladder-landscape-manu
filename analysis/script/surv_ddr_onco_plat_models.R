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
  broom::tidy(.)

dft_met_ddr_surv_sub <- dft_met_ddr_surv %>%
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
    race_ethnicity,
    female
  )

# Some problems with this as we have ECOG data that's about half missing.
coxph(
  Surv(
    time = fmr_fcpt_yrs,
    time2 = tt_os_first_met_reg_yrs,
    event = os_first_met_reg_status
  ) ~ ddr_before_entry + de_novo_met + md_ecog_imp_num + institution + race_ethnicity + female,
  data = dft_met_ddr_surv
)


dft_met_ddr_surv_sub <- dft_met_ddr_surv %>%
  select(
    fmr_fcpt_yrs,
    tt_os_first_met_reg_yrs,
    os_first_met_reg_status,
    ddr_before_entry,
    de_novo_met, 
    md_ecog_imp_num,
    institution,
    race_ethnicity,
    female
  ) %>%
  fastDummies::dummy_cols(
    select_columns = c("institution", "race_ethnicity"),
    remove_selected_columns = T, 
    remove_most_frequent_dummy = T
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

