library(purrr)
library(here)
library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source) # also load

dft_met_plat <- readr::read_rds(
  here('data', 'survival', 'first_line_platinum', 'data_first_line_plat.rds')
)

dir_out <- here('data', 'survival', 'first_line_platinum')


# Grab the stuff we want in the model, create indicator vars:
dft_met_plat %<>%
  select(
    # put the survival vars up front:
    reg_start_cpt_yrs,
    tt_os_g_yrs,
    os_g_status,
    # limit to the other variables we want for our model:
    carboplatin,
    bin_prev_plat,
    bin_prev_nonplat,
    md_ecog_imp_num,
    age_reg_start = dob_reg_start_yrs,
    de_novo_met,
    race_eth,
    female,
    institution
  ) %>%
  fastDummies::dummy_cols(
    select_columns = c("institution", "race_eth"),
    remove_selected_columns = T,
    remove_most_frequent_dummy = T
  )


surv_obj_os <- with(
  dft_met_plat,
  Surv(
    time = reg_start_cpt_yrs,
    time2 = tt_os_g_yrs,
    event = os_g_status
  )
)

# Some survival code I should probaly move later on:
dft_cox_univariate <- coxph(
  Surv(
    time = reg_start_cpt_yrs,
    time2 = tt_os_g_yrs,
    event = os_g_status
  ) ~
    carboplatin,
  data = dft_met_plat
) %>%
  broom::tidy(., conf.int = T)
# This gets save below in a big object with the other models.

# Run the complete case model
dft_cox_complete_cases <- coxph(
  Surv(
    time = reg_start_cpt_yrs,
    time2 = tt_os_g_yrs,
    event = os_g_status
  ) ~
    .,
  data = dft_met_plat
) %>%
  broom::tidy(., conf.int = T)
# This gets saved below along with the other models.

# We assume only md_ecog_imp_num to be missing.  Check with this:
# md.pattern(dft_met_plat, rotate.names = T)

# Do the multiple imputation, using every other covariate to predict the
#   missing ecog scores.

blocks <- construct.blocks(
  formulas = list(md_ecog_imp_num ~ .)
)
pm <- make.predictorMatrix(
  data = dplyr::select(
    dft_met_plat,
    -c(
      # just taking away the survival variables for predicting md_ecog_imp_sum
      reg_start_cpt_yrs,
      tt_os_g_yrs,
      os_g_status
    )
  ),
  blocks = blocks
)

# If you treat ECOG scores as an ordinal variable, which is probably
#   more correct, you get 'polr' as the default method.  I this seems
#   fine for imputation, but it doesn't make a ton of sense to me for
#   the model (3 terms).

mids_surv <- mice(
  data = dft_met_plat,
  maxit = 5,
  m = 15,
  seed = 234123498,
  predictorMatrix = pm,
  blocks = blocks,
  print = FALSE,
  # makes sense to me, gives similar results to pmm (default)
  method = 'midastouch'
)

gg_imp <- plot_gg_strip(mids_surv, var = "md_ecog_imp_num", pt_size = 1) +
  labs(y = "Last observed ECOG")

readr::write_rds(
  x = gg_imp,
  file = here(dir_out, 'gg_cox_imputation_md_ecog.rds')
)

# A bit painful here:  You can't use the data argument, or the "~." syntax.
#   If you do, then it will use the actual data without imputations.
#   I think there would be a way to autogen the formula and feed it in,
#   but it's only ~10 names so I'm going to type them out.
mira_surv <- with(
  mids_surv,
  coxph(
    Surv(
      time = reg_start_cpt_yrs,
      time2 = tt_os_g_yrs,
      event = os_g_status
    ) ~
      carboplatin +
        bin_prev_plat +
        bin_prev_nonplat +
        md_ecog_imp_num +
        age_reg_start +
        de_novo_met +
        female +
        race_eth_Asian +
        race_eth_Black +
        race_eth_race_other_unk +
        institution_DFCI +
        institution_UHN +
        institution_VICC
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

n_uni <- dft_met_plat %>%
  filter(!is.na(carboplatin)) %>%
  nrow(.)
n_cc <- dft_met_plat %>%
  filter(complete.cases(.)) %>%
  nrow(.)
n_mi <- dft_met_plat %>%
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


gg_cox_mod_compare <- ggplot(
  dat = mutate(dft_cox_all_mods, term = fct_rev(term)),
  aes(
    x = estimate,
    xmin = conf.low,
    xmax = conf.high,
    y = term,
    color = model_disp
  )
) +
  geom_vline(color = 'gray70', linewidth = 2, alpha = 0.5, xintercept = 0) +
  geom_pointrange(position = position_dodge2(width = 0.5), shape = 124) +
  theme_bw() +
  scale_color_highcontrast() +
  labs(y = NULL, x = "Cumulative log hazard ratio (95% CI)") +
  guides(color = guide_legend(title = NULL)) +
  theme(legend.position = "bottom")

readr::write_rds(
  gg_cox_mod_compare,
  here(dir_out, "gg_cox_mod_compare.rds")
)
