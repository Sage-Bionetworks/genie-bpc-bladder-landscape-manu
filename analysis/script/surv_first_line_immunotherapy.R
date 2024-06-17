
library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source) # also load
dft_ca_ind <- readr::read_rds(here('data', 'cohort', "ca_ind.rds"))
dft_reg <- readr::read_rds(here('data', 'cohort', "reg.rds"))
dft_cpt <- readr::read_rds(here('data', 'cohort', "cpt_aug.rds"))
dft_lot <- readr::read_rds(here('data', 'dmet', 'lines_of_therapy', 'lot.rds'))
dft_alt <- readr::read_rds(here('data', 'genomic','alterations.rds'))


dir_out <- here('data', 'survival', 'first_line_immuno')
fs::dir_create(dir_out)



lot_regex <- paste(
  c(
    'atezolizumab',
    'durvalumab',
    'cemiplimab',
    'nivolumab',
    'pembrolizumab',
    'ipilimumab',
    'tremelimumab'
  ),
  collapse = "|"
)

fl_immuno_regimens <- dft_lot %>% 
  filter(line_therapy %in% 1) %>% 
  count(regimen_drugs) %>%
  filter(str_detect(tolower(regimen_drugs), lot_regex)) 

readr::write_rds(
  fl_immuno_regimens,
  here(dir_out, "qualifying_regimen_count.rds")
)


dft_fl_io <- dft_lot %>%
  filter(line_therapy %in% 1) %>%
  mutate(
    immuno = factor(
      case_when(
        regimen_drugs %in% fl_immuno_regimens$regimen_drugs ~ "IO regimen",
        T ~ "No IO regimen"
      )
    )
  )

dft_fl_io <- dft_reg %>%
  select(
    record_id, ca_seq, regimen_number, 
    dx_reg_start_int_yrs,
    dx_reg_end_all_int_yrs,
    os_g_status, tt_os_g_yrs
  ) %>%
  left_join(
    dft_fl_io,
    .,
    by = c("record_id", "ca_seq", "regimen_number")
  )


dft_first_cpt <- get_first_cpt(
  dft_ca_ind,
  dft_cpt,
  include_sample_id = T
)

dft_fl_io <- left_join(
  dft_fl_io,
  dft_first_cpt,
  by = c('record_id', 'ca_seq')
) %>%
  mutate(
    reg_fcpt_yrs = dx_cpt_rep_yrs - dx_reg_start_int_yrs
  )



dft_fl_io %<>% 
  remove_trunc_gte_event(
    trunc_var = 'reg_fcpt_yrs',
    event_var = 'tt_os_g_yrs'
  )

dft_fl_io %<>% 
  mutate(reg_fcpt_yrs = ifelse(reg_fcpt_yrs < 0, 0, reg_fcpt_yrs))

surv_obj_os_fl_io <- with(
  dft_fl_io,
  Surv(
    time = reg_fcpt_yrs,
    time2 = tt_os_g_yrs,
    event = os_g_status
  )
)

gg_os_fl_io <- plot_one_survfit(
  dat = dft_fl_io,
  surv_form = surv_obj_os_fl_io ~ immuno,
  plot_title = "OS from first line therapy (metastatic)",
  plot_subtitle = "Adjusted for (independent) delayed entry"
)

readr::write_rds(
  dft_fl_io,
  here(dir_out, 'surv_first_line_immuno.rds')
)

readr::write_rds(
  gg_os_fl_io,
  here(dir_out, 'gg_os_fl_io.rds')
)
