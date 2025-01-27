# Creates an output versions of figures for GU ASCO 2025 submission poster.

library(purrr); library(fs); library(here);
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

dir_out <- here('output', 'fig', 'asco_gu_2025')

first_line_by_ddr <- readr::read_rds(
  file = here('data', 'survival', 'ddr_onco', 
              'gg_met_ddr_manu_no_rt.rds')
)

first_line_by_ddr_plat <- readr::read_rds(
  file = here('data', 'survival', 'ddr_onco', 
              'gg_met_ddr_plat_split_no_rt.rds')
)

second_line_by_ddr_io <- readr::read_rds(
  file = here('data', 'survival', 'ddr_onco', 
              'gg_os_2l_io_ddr_no_rt.rds')
)



gg_helper <- function(p, n) {
  ggsave(
    plot = p,
    filename = here(dir_out, paste0(n, '.pdf')),
    height = 2.5, width = 5, scale = 1.8
  )
  
  ggsave(
    plot = p,
    filename = here(dir_out, paste0(n, '.png')),
    height = 2.5, width = 5, scale = 1.8
  )
}

gg_helper(first_line_by_ddr, "1L_plat_by_ddr")
gg_helper(first_line_by_ddr_plat, "1L_plat_by_ddr_and_plat_base")
gg_helper(second_line_by_ddr_io, "2L_plat_by_ddr_and_plat_base")
