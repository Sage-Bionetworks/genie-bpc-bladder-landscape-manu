library(purrr)
library(here)
library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source) # also load

dft_met_ddr_surv <- readr::read_rds(
  here('data', 'survival', 'ddr_onco', 'met_ddr_surv.rds')
)

dir_out <- here('data', 'survival', 'ddr_onco')


# Just for clarity:
dft_met_ddr_surv %<>% rename(ddr_onco_alt = ddr_before_entry)


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
  ) ~
    ddr_onco_alt,
  data = dft_met_ddr_surv
) %>%
  broom::tidy(., conf.int = T)
# This gets save below in a big object with the other models.

dft_met_ddr_surv %<>%
  mutate(
    de_novo_met = if_else(
      dx_dmet_yrs > 0,
      1,
      0
    )
  )


dft_met_ddr_surv_mod_ready <- dft_met_ddr_surv %>%
  mutate(dob_reg_start_yrs = dob_reg_start_days / 365.25) %>%
  select(
    fmr_fcpt_yrs,
    tt_os_first_met_reg_yrs,
    os_first_met_reg_status,
    ddr_onco_alt,
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


# Some problems with this as we have ECOG data that's about half missing.
dft_cox_complete_cases <- coxph(
  Surv(
    time = fmr_fcpt_yrs,
    time2 = tt_os_first_met_reg_yrs,
    event = os_first_met_reg_status
  ) ~
    .,
  data = dft_met_ddr_surv_mod_ready
) %>%
  broom::tidy(., conf.int = T)
# This gets saved below along with the other models.

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
    )
  ),
  blocks = blocks
)

# If you treat ECOG scores as an ordinal variable, which is probably
#   more correct, you get 'polr' as the default method.  I this seems
#   fine for imputation, but it doesn't make a ton of sense to me for
#   the model (3 terms).

mids_surv <- mice(
  data = dft_met_ddr_surv_mod_ready,
  maxit = 5,
  m = 15,
  seed = 2341,
  predictorMatrix = pm,
  blocks = blocks,
  print = FALSE,
  # tended to be pretty similar to predictive mean matching:
  # method = 'cart'
  method = 'midastouch' # makes sense to me, gives similar results.
)

gg_imp <- plot_gg_strip(mids_surv, var = "md_ecog_imp_num", pt_size = 1) +
  labs(y = "Last observed ECOG")

readr::write_rds(
  x = gg_imp,
  file = here(dir_out, 'gg_cox_imputation_md_ecog.rds')
)

# # Put the multiply imputed datasets into a list column
# dat_res <- tibble(
#   dat = mice::complete(mids_surv, "all")
# )

# A bit painful here:  You can't use the data argument, or the "~." syntax.
#   If you do, then it will use the actual data without imputations.
#   I think there would be a way to autogen the formula and feed it in,
#   but it's only ~10 names so I'm going to type them out.
mira_surv <- with(
  mids_surv,
  coxph(
    Surv(
      time = fmr_fcpt_yrs,
      time2 = tt_os_first_met_reg_yrs,
      event = os_first_met_reg_status
    ) ~
      ddr_onco_alt +
        de_novo_met +
        age_reg_start +
        md_ecog_imp_num +
        female +
        institution_MSK +
        institution_UHN +
        institution_VICC +
        race_eth_Asian +
        race_eth_Black +
        race_eth_race_other_unk
  )
)

readr::write_rds(
  x = mira_surv,
  file = here(dir_out, 'cox_mult_imp_models.rds')
)


dft_cox_mult_imp <- pool(mira_surv) %>%
  summary(conf.int = T) %>%
  select(term, estimate, conf.low = `2.5 %`, conf.high = `97.5 %`, p.value) %>%
  mutate(level = 0.95) #alpha level, since I blurred that out above.

dft_cox_all_mods <- bind_rows(
  mutate(dft_cox_univariate, model = "uni"),
  mutate(dft_cox_complete_cases, model = "cc"),
  mutate(dft_cox_mult_imp, model = "mi")
) %>%
  mutate(
    term = str_replace_all(term, "TRUE$", ""),
    term = fct_inorder(term)
  )

readr::write_rds(
  dft_cox_all_mods,
  here(dir_out, "cox_tidy_all_models.rds")
)

n_uni <- dft_met_ddr_surv_mod_ready %>%
  filter(!is.na(ddr_onco_alt)) %>%
  nrow(.)
n_cc <- dft_met_ddr_surv_mod_ready %>%
  filter(complete.cases(.)) %>%
  nrow(.)
n_mi <- dft_met_ddr_surv_mod_ready %>%
  nrow(.)

dft_cox_all_mods %<>%
  mutate(
    model_disp = fct_inorder(
      case_when(
        model %in% "uni" ~ glue("Univariate (n={n_uni})"),
        model %in% "cc" ~ glue("Complete cases (n={n_cc})"),
        model %in% "mi" ~ glue("Multiple Imputation (n={n_mi})")
      )
    )
  )

# Adapt this code to use the new

# gg_cox_mod_compare <- ggplot(
#   dat = mutate(dft_cox_all_mods, term = fct_rev(term)),
#   aes(x = estimate, xmin = conf.low, xmax = conf.high, y = term,
#       color = model_disp)
# ) +
#   geom_vline(color = 'gray70', linewidth = 2, alpha = 0.5, xintercept = 0) +
#   geom_pointrange(position = position_dodge2(width = 0.5),
#                   shape = 124) +
#   theme_bw() +
#   scale_color_highcontrast() +
#   labs(y = NULL,
#        x = "Cumulative log hazard ratio (95% CI)") +
#   guides(color = guide_legend(title = NULL)) +
#   theme(legend.position = "bottom")

readr::write_rds(
  gg_cox_mod_compare,
  here(dir_out, "gg_cox_mod_compare.rds")
)


dft_alt_onco_ddr <- readr::read_rds(
  here('data', 'survival', 'ddr_onco', 'alt_onco_ddr.rds')
)

ercc2_altered <- dft_alt_onco_ddr %>%
  filter(hugo %in% 'ERCC2') %>%
  pull(record_id)

dft_met_ddr_surv %<>%
  mutate(
    ercc2 = record_id %in% ercc2_altered
  )

surv_ddr <- with(
  dft_met_ddr_surv,
  Surv(
    time = fmr_fcpt_yrs,
    time2 = tt_os_first_met_reg_yrs,
    event = os_first_met_reg_status
  )
)

survfit(surv_ddr ~ ercc2, data = dft_met_ddr_surv)
coxph(surv_ddr ~ ercc2, data = dft_met_ddr_surv)
