# _interact in the title refers to the goal of looking at the interaction between
#   DDR and platinum therapy.  In other words, do DDR people have a better 
#   response to platinum, or do they just happen to do better on first line therapy?

library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

dir_out <- here('data', 'survival', 'ddr_onco_1L')

met_ddr_surv <- readr::read_rds(
  here('data', 'survival', 'ddr_onco_1L', 'met_ddr_surv_all_1L.rds')
)

met_ddr_surv %<>% rename(ddr_onco_alt = ddr_before_entry)

# For this analysis we'll group regimens a bit more:
met_ddr_surv %<>%
  mutate(
    regimen_cat_km = forcats::fct_collapse(
      regimen_cat,
      Other = c("Taxane monotherapy",
                "Pemetrexed monotherapy",
                "Other")
    )
  )

# Interval from regimen to first cancer panel test - make zero if neg.
met_ddr_surv %<>% mutate(reg_fcpt_yrs = ifelse(reg_fcpt_yrs < 0, 0, reg_fcpt_yrs))
# These people wouldn't exist if we had proper cohort entry, but we can't do much with them:
met_ddr_surv %<>% 
  remove_trunc_gte_event(
    trunc_var = 'reg_fcpt_yrs',
    event_var = 'tt_os_g_yrs'
  )

met_ddr_surv %<>%
  mutate(
    ddr_lab = factor(
      case_when(ddr_onco_alt ~ "DDR+",
                T ~ "DDR-")
    )
  )
                

model_bundle <- list(data_km = met_ddr_surv)


individual_km <- met_ddr_surv %>%
  nest(.by = c(regimen_cat_km, ddr_lab)) %>%
  mutate(
    km_tidy = purrr::map(
      .x = data,
      .f = \(x) km_tidy_no_covariate(
        dat = x,
        trunc_time = 'reg_fcpt_yrs',
        event_time = 'tt_os_g_yrs',
        event_ind = 'os_g_status'
      )
    )
  )

individual_km %<>% 
  select(-data) %>%
  unnest(km_tidy)

individual_km %<>%
  group_by(regimen_cat_km) %>%
  mutate(drug_n = sum(records)) %>%
  ungroup() %>%
  arrange(desc(drug_n), regimen_cat_km, desc(ddr_lab)) %>%
  mutate(ddr_drug_axis = fct_inorder(
    glue('(n={records}) {regimen_cat_km}; {ddr_lab}')
  ))

model_bundle <- c(model_bundle, list(individual_km = individual_km))


gg_km_1L <- plot_km_forest(
  individual_km,
  y = "ddr_drug_axis",
  plot_infinite = T
) +
  labs(
    title = "KM median/CI for 1L",
    subtitle = "Each estimate is specific to one drug/ddr group",
    y = NULL,
    x = "Median overall survival (years)"
  )

model_bundle <- c(model_bundle, list(gg_km_1L = gg_km_1L))







# Cox model fitting
met_ddr_surv_mod_ready <- met_ddr_surv %>%
  mutate(
    platinum = regimen_cat_km %in% c("Carboplatin-based", "Cisplatin-based") 
  ) %>%
  select(
    reg_fcpt_yrs,
    tt_os_g_yrs,
    os_g_status,
    platinum,
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


model_bundle <- c(model_bundle, data_lm = list(met_ddr_surv_mod_ready))




# Do the multiple imputation, using every other covariate to predict the
#   missing ecog scores.

blocks <- construct.blocks(
  formulas = list(md_ecog_imp_num ~ .)
)
pm <- make.predictorMatrix(
  data = dplyr::select(
    met_ddr_surv_mod_ready, 
    -c(
      # just taking away the survival variables for predicting md_ecog_imp_sum
      reg_fcpt_yrs,
      tt_os_g_yrs,
      os_g_status
    )),
  blocks = blocks
)

# If you treat ECOG scores as an ordinal variable, which is probably
#   more correct, you get 'polr' as the default method.  I this seems
#   fine for imputation, but it doesn't make a ton of sense to me for
#   the model (3 terms).

mids_surv <- mice(
  data = met_ddr_surv_mod_ready, 
  maxit = 5, m = 15, seed = 2341, 
  predictorMatrix = pm, 
  blocks = blocks,
  print = FALSE,
  # tended to be pretty similar to predictive mean matching:
  # method = 'cart'
  method = 'midastouch' # makes sense to me, gives similar results.
)

gg_imp <- plot_gg_strip(
  mids_surv, var = "md_ecog_imp_num",
  pt_size = 1
) + 
  labs(y = "Last observed ECOG")

model_bundle <- c(model_bundle, gg_impute = list(gg_imp))

# A bit painful here:  You can't use the data argument, or the "~." syntax.
#   If you do, then it will use the actual data without imputations.
#   I think there would be a way to autogen the formula and feed it in,
#   but it's only ~10 names so I'm going to type them out.
mira_surv <- with(
  mids_surv,
  coxph(
    Surv(
      time = reg_fcpt_yrs,
      time2 = tt_os_g_yrs,
      event = os_g_status
    ) ~ ddr_onco_alt * platinum +
      de_novo_met + 
      age_reg_start + 
      md_ecog_imp_num + 
      female + 
      institution_DFCI + 
      institution_UHN + 
      institution_VICC + 
      race_eth_Asian + 
      race_eth_Black + 
      race_eth_race_other_unk
  )
)

model_bundle <- c(model_bundle, fit_imputation = list(mira_surv))

cox_mult_imp <- pool(mira_surv) %>% 
  summary(conf.int = T) %>% 
  select(term, estimate, conf.low = `2.5 %`, conf.high = `97.5 %`, p.value) %>%
  mutate(level = 0.95) #alpha level, since I blurred that out above.

model_bundle <- c(model_bundle, tidy_imputed = list(cox_mult_imp))

gg_cox_mod_imp <- ggplot(
  dat = mutate(cox_mult_imp, term = fct_rev(term)),
  aes(x = estimate, xmin = conf.low, xmax = conf.high, y = term)
) + 
  geom_vline(color = 'gray70', linewidth = 2, alpha = 0.5, xintercept = 0) + 
  geom_pointrange(position = position_dodge2(width = 0.5),
                  shape = 124) + 
  theme_bw() + 
  scale_color_highcontrast() +
  labs(
    title = "1L model with Plat:DDR interaction",
    subtitle = "Multiple imputation used",
    y = NULL, 
    x = "Cumulative log hazard ratio (95% CI)") + 
  guides(color = guide_legend(title = NULL)) +
  theme(legend.position = "bottom")

model_bundle <- c(model_bundle, gg_imputed = list(gg_cox_mod_imp))


readr::write_rds(
  model_bundle,
  here(dir_out, 'model_interact_bundle.rds')
)


