
library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source) # also load


dir_out <- here('data', 'survival', 'ddr_onco')

dft_met_ddr_surv <- read_rds(here(dir_out, 'met_ddr_surv.rds'))
  
# Create the basic survival plots

dft_met_ddr_surv %<>% 
  remove_trunc_gte_event(
    trunc_var = 'fmr_fcpt_yrs',
    event_var = 'tt_os_first_met_reg_yrs'
  )

dft_met_ddr_surv %<>% mutate(fmr_fcpt_yrs = ifelse(fmr_fcpt_yrs < 0, 0, fmr_fcpt_yrs))

surv_obj_os_fmr <- with(
  dft_met_ddr_surv,
  Surv(
    time = fmr_fcpt_yrs,
    time2 = tt_os_first_met_reg_yrs,
    event = os_first_met_reg_status
  )
)

dft_met_ddr_surv %<>%
  mutate(
    ddr_disp = case_when(
      ddr_before_entry ~ "Oncogenic DDR",
      T ~ "No Onco. DDR"
    )
  )

gg_os_fmr_ddr <- plot_one_survfit(
  dat = dft_met_ddr_surv,
  surv_form = surv_obj_os_fmr ~ ddr_disp,
  plot_title = "OS from first line platinum chemo",
  plot_subtitle = "Adjusted for (independent) delayed entry"
)

readr::write_rds(
  gg_os_fmr_ddr,
  file = here(dir_out, "gg_met_ddr.rds")
)


# Make a truncated version:

gg_os_fmr_ddr_manu <- plot_one_survfit(
  dat = dft_met_ddr_surv,
  surv_form = surv_obj_os_fmr ~ ddr_disp,
  plot_title = "OS from first line platinum chemo",
  x_breaks = seq(0, 100, by = 0.5)
) + 
  coord_cartesian(xlim = c(0,5)) +
  theme(plot.title.position = 'panel')

gg_os_fmr_ddr_manu_no_rt <- plot_one_survfit_no_risktable(
  dat = dft_met_ddr_surv,
  surv_form = surv_obj_os_fmr ~ ddr_disp,
  plot_title = "OS from first line platinum chemo",
  x_breaks = seq(0, 100, by = 0.5)
) + 
  coord_cartesian(xlim = c(0,5)) +
  theme(plot.title.position = 'panel')

readr::write_rds(
  gg_os_fmr_ddr_manu,
  file = here(dir_out, "gg_met_ddr_manu.rds")
)

readr::write_rds(
  gg_os_fmr_ddr_manu_no_rt,
  file = here(dir_out, "gg_met_ddr_manu_no_rt.rds")
)








# Request:  Also do this split by carboplatin and cisplatin.
dft_met_ddr_surv %<>%
  mutate(
    carbo = str_detect(regimen_drugs, "arboplatin"),
    cis = str_detect(regimen_drugs, 'isplatin')
  ) %>%
  filter(!(carbo & cis)) %>% # one case, not real (0 length)
  mutate(
    ddr_plat_comb = case_when(
      ddr_before_entry & carbo ~ "Carbo, DDR+",
      ddr_before_entry & cis ~ "Cis, DDR+",
      carbo ~ "Carbo, DDR-",
      cis ~ "Cis, DDR-"
    )
  ) %>%
  mutate(ddr_plat_comb = factor(ddr_plat_comb))

surv_obj_os_fmr <- with(
  dft_met_ddr_surv,
  Surv(
    time = fmr_fcpt_yrs,
    time2 = tt_os_first_met_reg_yrs,
    event = os_first_met_reg_status
  )
)

gg_os_fmr_ddr_manu_plat_split <- plot_one_survfit(
  dat = dft_met_ddr_surv,
  surv_form = surv_obj_os_fmr ~ ddr_plat_comb,
  plot_title = "OS from first line platinum chemo",
  x_breaks = seq(0, 100, by = 0.5),
  pal = c("#ee99aa", "#994455", "#6699cc", "#004488")
) + 
  coord_cartesian(xlim = c(0,5)) +
  theme(plot.title.position = 'panel')

gg_os_fmr_ddr_manu_plat_split_no_rt <- plot_one_survfit_no_risktable(
  dat = dft_met_ddr_surv,
  surv_form = surv_obj_os_fmr ~ ddr_plat_comb,
  plot_title = "OS from first line platinum chemo",
  x_breaks = seq(0, 100, by = 0.5),
  pal = c("#ee99aa", "#994455", "#6699cc", "#004488")
) + 
  coord_cartesian(xlim = c(0,5)) +
  theme(plot.title.position = 'panel')

gg_os_fmr_ddr_manu_plat_split_no_rt


readr::write_rds(
  gg_os_fmr_ddr_manu_plat_split,
  file = here(dir_out, "gg_met_ddr_plat_split.rds")
)

readr::write_rds(
  gg_os_fmr_ddr_manu_plat_split_no_rt,
  file = here(dir_out, "gg_met_ddr_plat_split_no_rt.rds")
)




