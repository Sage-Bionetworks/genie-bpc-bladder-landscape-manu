# Yet another iteration of this question, rephrased by Michal and Neil for the
# 2024 ASCO urinary abstract deadline.

library(purrr)
library(here)
library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source) # also load

met_plat_1l <- readr::read_rds(
  here('data', 'survival', 'ddr_onco_1L', 'met_ddr_surv_plat_1L.rds')
)

met_plat_1l %<>%
  remove_trunc_gte_event(
    trunc_var = 'reg_fcpt_yrs',
    event_var = 'tt_os_g_yrs'
  )

met_plat_1l %<>%
  mutate(reg_fcpt_yrs = ifelse(reg_fcpt_yrs < 0, 0, reg_fcpt_yrs))


met_plat_1l_grouped <- bind_rows(
  (met_plat_1l %>% mutate(analysis_group = "any_plat")),
  (met_plat_1l %>%
    filter(str_detect(regimen_drugs, "Cisplatin")) %>%
    mutate(analysis_group = "cisplatin_based")),
  (met_plat_1l %>%
    filter(str_detect(regimen_drugs, "Carboplatin")) %>%
    mutate(analysis_group = "carboplatin_based"))
)

met_plat_1l_grouped %<>%
  nest(.by = analysis_group)


met_plat_1l_grouped %<>%
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
              time2 = tt_os_g_yrs,
              event = os_g_status
            )
          ) ~
            ddr_before_entry
        ) %>%
          broom::tidy(., conf.int = T, exponentiate = F) # will exp later.
      }
    )
  )


met_plat_1l_grouped %<>%
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
              time2 = tt_os_g_yrs,
              event = os_g_status
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

# Removed a printout here.

met_plat_1l_grouped %>%
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

# This is the one that ended up in Summer summit 2025 presentation:
met_plat_1l_grouped %>%
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
