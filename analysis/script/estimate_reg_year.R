
library(fs); library(purrr); library(here);
purrr::walk(.x = fs::dir_ls('R'), .f = source)

read_wrap <- function(p) {
  read_rds(
    file = here("data", 'cohort', paste0(p, '.rds'))
  )
}

dft_pt <- read_wrap("pt")
dft_ca_ind <- read_wrap("ca_ind")
dft_reg <- read_wrap("reg")

dft_reg_yr <- dft_reg %>%
  select(record_id, ca_seq, regimen_number, dx_reg_start_int_yrs) %>%
  # add in the data we need for birthdays.
  left_join(
    .,
    select(dft_pt, record_id, birth_year),
    by = "record_id"
  ) %>%
  left_join(
    .,
    select(dft_ca_ind, record_id, dob_ca_dx_yrs),
    by = "record_id"
  )


dft_reg_yr %<>%
  mutate(
    reg_start_date_jan1_bd = ymd(paste0(birth_year,"-01-01")) +
      dyears(dob_ca_dx_yrs) +
      dyears(dx_reg_start_int_yrs),
    reg_start_date_dec31_bd = ymd(paste0(birth_year,"-12-31")) +
      dyears(dob_ca_dx_yrs) +
      dyears(dx_reg_start_int_yrs),
  ) %>%
  mutate(
    reg_start_yr_jan1_bd = year(reg_start_date_jan1_bd),
    reg_start_yr_dec31_bd = year(reg_start_date_dec31_bd)
  )

readr::write_rds(
  x = dft_reg_yr,
  file = here('data', 'cohort', 'reg_yr_case_level.rds')
)

one_iter_yr <- function(s, dat) {
  set.seed(s)
  dat %>%
    mutate(
      add_days = sample.int(365, n(), replace = T) - 1,
      event_date = reg_start_date_jan1_bd + ddays(add_days),
      event_yr = year(event_date)
    ) %>%
    count(event_yr)
}

many_iter_yr <- function(input_dat, n_rep) {
  rtn <- tibble(id = 1:n_rep) %>%
    mutate(s = sample.int(n_rep*10^3, size = n(), replace = T)) %>%
    mutate(dat = purrr::map(.x = s,
                            .f = one_iter_yr,
                            dat = input_dat)) %>%
    unnest(dat) %>%
    mutate(event_yr = factor(event_yr)) %>%
    tidyr::complete(id, event_yr, fill = list(n = 0))
}

sum_reps <- function(rep_dat) {
  rep_dat %>%
    group_by(event_yr) %>%
    summarize(
      median_n = median(n, na.rm = T),
      lcb_n = quantile(n, probs = c(0.025)),
      ucb_n = quantile(n, probs = c(0.975)),
      # yr_ct_sanity = n(), # does not not need to be output.
      .groups = "drop"
    )
}

reps <- many_iter_yr(input_dat = dft_reg_yr, n_rep = 10^3)
rep_sum <- sum_reps(reps)

readr::write_rds(
  x = rep_sum,
  file = here('data', 'cohort', 'reg_yr_estimates.rds')
)




# Update:  got a request to do this again for met regimens only.
dft_reg_met <- readr::read_rds(here('data', 'dmet', 'reg_start_gte_dmet.rds'))

dft_reg_yr_met <- dft_reg_met %>% 
  select(record_id, ca_seq, regimen_number) %>%
  distinct(.) %>%
  left_join(
    .,
    dft_reg_yr,
    by = c("record_id", "ca_seq", "regimen_number"),
    relationship = "one-to-many"
  ) 

dft_reg_yr_non_met <- anti_join(
  dft_reg_yr, 
  dft_reg_yr_met,
  by = c("record_id", "ca_seq", "regimen_number")
)

reps_met <- many_iter_yr(input_dat = dft_reg_yr_met, n_rep = 10^3)
rep_sum_met <- sum_reps(reps_met)

reps_non_met <- many_iter_yr(input_dat = dft_reg_yr_non_met, n_rep = 10^3)
rep_sum_non_met <- sum_reps(reps_non_met)

dft_rep_comb <- bind_rows(
  mutate(rep_sum, type = "All"),
  mutate(rep_sum_met, type = "Metastatic"),
  mutate(rep_sum_non_met, type = "Not metastatic")
) %>%
  mutate(type = fct_inorder(type))
  
readr::write_rds(
  x = dft_rep_comb,
  file = here('data', 'cohort', 'reg_yr_estimates_three_ways.rds')
)

  
  
