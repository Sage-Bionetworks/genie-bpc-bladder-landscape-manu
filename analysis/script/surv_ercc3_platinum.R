# Create analysis dataset for a specific case identified by our physicians:
# Platinum based chemotherapy with gemcitabine, does it increase survival to 
#   have ERCC3 alterations?

library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source) 

dft_ca_ind <- readr::read_rds(here('data', 'cohort', "ca_ind.rds"))
dft_reg <- readr::read_rds(here('data', 'cohort', "reg.rds"))
dft_cpt <- readr::read_rds(here('data', 'cohort', "cpt_aug.rds"))
dft_alt <- readr::read_rds(here('data', 'genomic','alterations.rds'))

dir_out <- here('data', 'survival', 'ercc3_plat')
fs::dir_create(dir_out)



# The vector of regimens we'll consider to be platinum based chemo.
vec_gem_plat <- c(
  "Cisplatin, Gemcitabine Hydrochloride",
  "Carboplatin, Gemcitabine Hydrochloride"
)

dft_ercc3_alt <- dft_alt %>%
  filter(hugo %in% "ERCC3")

dft_ercc3_alt <- dft_cpt %>% 
  select(
    cpt_genie_sample_id, record_id, ca_seq, 
    dx_path_proc_cpt_yrs, dx_cpt_rep_yrs
  ) %>%
  left_join(
    dft_ercc3_alt,
    .,
    by = c(sample_id = "cpt_genie_sample_id")
  )

readr::write_rds(
  x = dft_ercc3_alt,
  file = here(dir_out, 'ercc3_alt.rds')
)

dft_first_cpt <- get_first_cpt(
  dft_ca_ind,
  dft_cpt,
  include_sample_id = T
)

dft_ercc3_sum <- dft_ercc3_alt %>%
  filter(sample_id %in% dft_first_cpt$cpt_genie_sample_id) %>%
  group_by(record_id, ca_seq) %>%
  summarize(
    gene_ercc3 = T,
    .groups = "drop"
  ) 

# First platinum based regimen after at least one cancer panel test.
dft_plat_first <- dft_reg %>%
  filter(regimen_drugs %in% vec_gem_plat) %>%
  group_by(record_id, ca_seq) %>%
  arrange(regimen_number) %>%
  slice(1) %>%
  ungroup(.) %>%
  select(
    record_id, ca_seq, regimen_drugs, 
    dx_reg_start_int, dx_reg_end_all_int,
    os_g_status, tt_os_g_days
  )

dft_plat_first %<>%
  left_join(
    .,
    dft_ercc3_sum,
    by = c("record_id", "ca_seq")
  ) %>%
  replace_na(replace = list(gene_ercc3 = F))


dft_plat_first <- dft_first_cpt %>%
  left_join(
    select(
      dft_cpt, record_id, ca_seq, cpt_genie_sample_id, 
      dx_path_proc_cpt_days, dx_cpt_rep_days,
    ),
    .,
    by = c("record_id", "ca_seq", "cpt_genie_sample_id")
  ) %>%
  left_join(
    dft_plat_first,
    .,
    by = c("record_id", "ca_seq")
  )
    
dft_plat_ercc3 <- dft_plat_first       

dft_plat_ercc3 %<>% 
  mutate(
    reg_cpt_rep_days = dx_cpt_rep_days - dx_reg_start_int,
    reg_cpt_rep_yrs = reg_cpt_rep_days / 365.25,
    tt_os_g_yrs = tt_os_g_days / 365.25
  ) 

dft_plat_ercc3 %<>%
  mutate(
    ercc3_disp = factor(case_when(
      gene_ercc3 ~ "ERCC3 altered",
      T ~ "No ERCC3 alteration"
    ))
  )





dft_plat_ercc3 %<>% 
  remove_trunc_gte_event(
    trunc_var = 'reg_cpt_rep_yrs',
    event_var = 'tt_os_g_yrs'
  )

dft_plat_ercc3 %<>% 
  mutate(reg_cpt_rep_yrs = ifelse(reg_cpt_rep_yrs < 0, 0, reg_cpt_rep_yrs))

surv_obj_os_plat_ercc3 <- with(
  dft_plat_ercc3,
  Surv(
    time = reg_cpt_rep_yrs,
    time2 = tt_os_g_yrs,
    event = os_g_status
  )
)

gg_os_plat_ercc3 <- plot_one_survfit(
  dat = dft_plat_ercc3,
  surv_form = surv_obj_os_plat_ercc3 ~ ercc3_disp,
  plot_title = "OS from first CarboGem or GemCis regimen",
  plot_subtitle = "Adjusted for (independent) delayed entry"
)

readr::write_rds(
  dft_plat_ercc3,
  here(dir_out, 'surv_dat_plat_ercc3.rds')
)

readr::write_rds(
  gg_os_plat_ercc3,
  here(dir_out, 'gg_surv_ercc3.rds')
)
