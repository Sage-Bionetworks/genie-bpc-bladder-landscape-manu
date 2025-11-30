library(purrr)
library(here)
library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

manu_out_dir_fig3 <- here('output', 'manu', 'fig3')

# These are already created:
gg_fig_3a <- readr::read_rds(
  here('data', 'survival', 'ddr_onco_1L', 'gg_met_ddr_manu.rds')
)
gg_fig_3b <- readr::read_rds(
  here('data', 'survival', 'ddr_onco_1L', 'gg_met_ddr_plat_split.rds')
)

gg_fig_3a <- gg_fig_3a + labs(title = NULL)
gg_fig_3b <- gg_fig_3b + labs(title = NULL)

manu_plot_save_helper(
  dir = manu_out_dir_fig3,
  gg_fig_3a,
  name = 'fig_3a_all_plat',
  height = 6,
  width = 8
)

manu_plot_save_helper(
  dir = manu_out_dir_fig3,
  gg_fig_3b,
  name = 'fig_3b_cis_carbo_split',
  height = 6,
  width = 8,
)

# If you do want to combine them into one output figure remember that ggsurvfit_build() works well for this:
# cowplot::plot_grid(
#   ggsurvfit_build(gg_fig_3a),
#   ggsurvfit_build(gg_fig_3b),
#   ...
# )

cli_abort(
  "still need the table outputs for these figures which are pasted up top in the 'plan'"
)
