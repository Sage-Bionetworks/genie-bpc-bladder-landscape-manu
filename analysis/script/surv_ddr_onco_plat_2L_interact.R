# _interact in the title refers to the goal of looking at the interaction between
#   DDR and platinum therapy.  In other words, do DDR people have a better 
#   response to platinum, or do they just happen to do better on first line therapy?

library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

dir_out <- here('data', 'survival', 'ddr_onco_2L')

met_ddr_surv <- readr::read_rds(
  here('data', 'survival', 'ddr_onco_2L', 'met_ddr_surv_all_2L.rds')
)

met_ddr_surv %<>% rename(ddr_onco_alt = ddr_before_entry)

# For this analysis we'll group regimens a bit more:
met_ddr_surv %<>%
  mutate(
    regimen_cat_km = forcats::fct_collapse(
      regimen_cat,
      Other = c(
        "IO (other)",
        "Taxane monotherapy",
        "Pemetrexed monotherapy",
        "Other"),
      `Platinum-based` = c(
        "Cisplatin-based",
        "Carboplatin-based"
      )
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


model_bundle <- list(data = met_ddr_surv)


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


gg_km_2L <- plot_km_forest(
  individual_km,
  y = "ddr_drug_axis",
  plot_infinite = T
) +
  labs(
    title = "KM median/CI for 2L",
    subtitle = "Each estimate is specific to one drug/ddr group",
    y = NULL,
    x = "Median overall survival (years)"
  )

model_bundle <- c(model_bundle, list(gg_km_2L))

readr::write_rds(
  model_bundle,
  here(dir_out, 'model_interact_bundle.rds')
)

# First output: Kaplan Meier estimates for
