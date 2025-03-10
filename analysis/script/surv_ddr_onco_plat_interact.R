# _interact in the title refers to the goal of looking at the interaction between
#   DDR and platinum therapy.  In other words, do DDR people have a better 
#   response to platinum, or do they just happen to do better on first line therapy?

library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

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

model_bundle <- list(data = met_ddr_surv)


test <- met_ddr_surv %>%
  nest(.by = c(regimen_cat_km, ddr_onco_alt)) %>%
  slice(1) %>%
  pull(data) %>%
  `[[`(.,1)

get_no_covariate_km(dat = test,
                    trunc_time = 'reg_fcpt_yrs',
                    event_time = 'tt_os_g_yrs',
                    event_ind = 'os_g_status')

surv_obj <- with(
  test,
  Surv(
    time = reg_fcpt_yrs,
    time2 = tt_os_g_yrs,
    event = os_g_status
  )
)

survfit(surv_obj ~ 1) %>%
  broom::glance(.)

km_tidy_no_covariate <- function(
    dat,
    event_time,
    event_ind,
    trunc_time = NULL
) {
  if (is.null(trunc)) {
    surv_obj <- Surv(time = dat[[event_time]], event = dat[[event_ind]])
  } else {
    surv_obj <- with(dat, Surv(time = dat[[trunc_time]],
                               time2 = dat[[event_time]],
                               event = dat[[event_ind]]))
  }
  survfit(surv_obj ~ 1, data = dat) %>%
    broom::glance(.)
  
}
    
       
# First output: Kaplan Meier estimates for
