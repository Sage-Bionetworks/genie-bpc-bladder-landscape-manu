# Creates an output version of Figure 1 for the paper.

library(purrr); library(fs); library(here);
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

dir_out <- here('output', 'fig', 'manu')

os_by_stage <- readr::read_rds(
    file = here('data', 'survival', 'os_dx_by_stage_manu.rds')
)

first_line_platinum <- readr::read_rds(
  file = here('data', 'survival', 'first_line_platinum', 'gg_first_line_platinum.rds')
)

first_line_by_ddr <- readr::read_rds(
  file = here('data', 'survival', 'ddr_onco', 
              'gg_met_ddr_manu.rds')
)

first_line_by_ddr_plat <- readr::read_rds(
  file = here('data', 'survival', 'ddr_onco', 
              'gg_met_ddr_plat_split.rds')
)


fig1_comb <- plot_grid(
  ggsurvfit_build(os_by_stage),
  ggsurvfit_build(first_line_platinum),
  ggsurvfit_build(first_line_by_ddr),
  ggsurvfit_build(first_line_by_ddr_plat),
  labels = "AUTO",
  align = 'hv',
  axis = 't'
)



ggsave(
  plot = fig1_comb,
  filename = here(dir_out, 'fig1.pdf'),
  height = 5, width = 7, scale = 1.8
)

ggsave(
  plot = fig1_comb,
  filename = here(dir_out, 'fig1.png'),
  height = 5, width = 7, scale = 1.8
)