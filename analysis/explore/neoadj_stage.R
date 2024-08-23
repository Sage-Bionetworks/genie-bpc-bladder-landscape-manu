library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

dft_ca_ind <- readr::read_rds(here('data', 'cohort', "ca_ind.rds"))

ex <- dft_ca_ind %>% 
  filter(ca_tx_pre_path_stage %in% c("Yes", "No")) %>%
  filter(stage_dx_iv %in% "Stage IV") %>%
  filter(institution %in% "MSK") %>%
  group_by(best_ajcc_stage_cd, ca_tx_pre_path_stage) %>%
  slice(1:2) %>%
  ungroup(.)

ex %<>%
  select(record_id, ca_seq, best_ajcc_stage_cd, stage_dx_iv, ca_tx_pre_path_stage)

readr::write_csv(
  ex,
  file = here('analysis', 'explore', 'neoadj_stage_4.csv')
)
