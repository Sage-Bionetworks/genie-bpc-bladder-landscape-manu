#

library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source) # also loads lots of packages.

surv_desc_fp <- here('data', 'survival')

read_wrap_clin <- function(p) {
  read_rds(file = here("data", 'cohort', p))
}

dft_pt <- read_wrap_clin("pt.rds")
dft_ca_ind <- read_wrap_clin("ca_ind.rds")
# we're taking the augmented version, which has TMB columns added.  The existing info is the same.
dft_cpt <- read_wrap_clin("cpt_aug.rds")

dft_lot <- readr::read_rds(
  here('data', 'dmet', 'lines_of_therapy', 'lot.rds')
)

dft_cases <- dft_lot %>%
  filter(line_therapy %in% 1) %>%
  filter(
    # so not including things like GemCis + paclitaxel which seem common.
    regimen_drugs %in% c(
      "Cisplatin, Gemcitabine Hydrochloride",
      "Carboplatin, Gemcitabine Hydrochloride"
    )
  )

dft_cases <- left_join(
  dft_cases,
  select(dft_reg, record_id, ca_seq, regimen_number,
         dx_reg_start_int_yrs, tt_os_g_yrs, os_g_status),
  by = c("record_id", "ca_seq", "regimen_number"),
  relationship = "one-to-one"
)

dft_first_cpt <- get_first_cpt(
  ca_ind_dat = dft_ca_ind,
  cpt_dat = dft_cpt
)

dft_cases <- left_join(
  dft_cases,
  dft_first_cpt,
  by = c("record_id", "ca_seq")
)

dft_cases %<>%
  mutate(
    reg_start_cpt_yrs = dx_cpt_rep_yrs - dx_reg_start_int_yrs
  ) %>%
  select(-c(dx_cpt_rep_yrs, dx_reg_start_int_yrs)) # avoid confusion

dft_cases %<>% remove_trunc_gte_event(
  trunc_var = 'reg_start_cpt_yrs',
  event_var = 'tt_os_g_yrs'
)

surv_obj_cases <- with(
  dft_cases,
  Surv(
    time = reg_start_cpt_yrs,
    time2 = tt_os_g_yrs,
    event = os_g_status
  )
)

gg_first_line_platinum <- plot_one_survfit(
  dat = dft_cases,
  surv_form = surv_obj_cases ~ regimen_drugs,
  plot_title = "OS, from start of first line platinum-based chemotherapy",
  plot_subtitle = "Adjusted for truncation assuming independence",
  x_exp = 0.1,
  x_breaks = 0:100
) + 
  add_quantile(y_value = 0.5, linetype = 81, alpha = 0.75)

readr::write_rds(
  gg_first_line_platinum,
  here('data', 'survival', 'first_line_platinum', 'gg_first_line_platinum.rds')
)

