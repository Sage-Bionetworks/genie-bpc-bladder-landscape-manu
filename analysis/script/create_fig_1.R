# Creates an output version of Figure 1 for the paper.

library(purrr); library(fs); library(here);
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

os_by_stage <- readr::read_rds(
    file = here('data', 'survival', 'os_dx_by_stage_ss24.rds')
)


plot_grid(
  ggsurvfit_build(os_by_stage),
  
  labels = "AUTO"
)
