library(purrr)
library(here)
library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

manu_out_dir_fig3 <- here('output', 'manu', 'fig3')

# These are already created:
gg_fig_3a <- readr::read_rds(
  here('data', 'survival', 'ddr_onco_1L', 'gg_met_ddr_manu.rds')
)
gg_fig_3x <- readr::read_rds(
  here('data', 'survival', 'ddr_onco_1L', 'gg_met_ddr_manu_iyer.rds')
)

manu_plot_save_helper(
  dir = manu_out_dir_fig3,
  cowplot::plot_grid(
    ggsurvfit_build(gg_fig_3a),
    ggsurvfit_build(gg_fig_3x),
    ncol = 2
  ),
  name = 'fig_3x_all_plat_iyer',
  height = 6,
  width = 16,
  as_jpeg = FALSE,
  as_png = FALSE
)
