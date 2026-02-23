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

fig3a_table_inlay <- readr::read_rds(
  here('data', 'survival', 'ddr_onco_1L', '3a_tab_inlay.rds')
)

fig3b_table_inlay <- readr::read_rds(
  here('data', 'survival', 'ddr_onco_1L', '3b_tab_inlay.rds')
)

fig3a_inlay_raster <- fig3a_table_inlay %>%
  rename(` ` = left_head) %>%
  flextable(.) %>%
  fontsize(size = 8, part = "all") %>%
  padding(padding = 3, part = "all") %>%
  flextable::merge_at(i = 3, j = 2:3) %>%
  autofit(.) %>%
  theme_booktabs(.) %>%
  bg(bg = "white", part = 'all') %>%
  gen_grob(., fit = 'fixed', scaling = 'fixed')

gg_fig_3a_with_inlay <- gg_fig_3a +
  annotation_custom(
    fig3a_inlay_raster,
    xmin = 1.6,
    xmax = 5.5,
    ymin = 0.75,
    ymax = 0.98
  )

manu_plot_save_helper(
  dir = manu_out_dir_fig3,
  gg_fig_3a_with_inlay,
  name = 'fig_3a_all_plat_inlay_table',
  height = 6,
  width = 8
)


fig3b_inlay_raster <- fig3b_table_inlay %>%
  rename(`  ` = left_grp_head, ` ` = left_head) %>%
  flextable(.) %>%
  fontsize(size = 8, part = "all") %>%
  padding(padding = 3, part = "all") %>%
  flextable::merge_v(j = 1) %>%
  flextable::valign(valign = "top") %>%
  flextable::merge_at(i = 3, j = 3:4) %>%
  flextable::merge_at(i = 6, j = 3:4) %>%
  autofit(.) %>%
  theme_booktabs(.) %>%
  bg(bg = "white", part = "all") %>%
  gen_grob(., fit = 'fixed', scaling = 'fixed')

gg_fig_3b_with_inlay <- gg_fig_3b +
  annotation_custom(
    fig3b_inlay_raster,
    xmin = 1.6,
    xmax = 5.5,
    ymin = 0.75,
    ymax = 0.98
  )

manu_plot_save_helper(
  dir = manu_out_dir_fig3,
  gg_fig_3b_with_inlay,
  name = 'fig_3b_cis_carbo_split_inlay_table',
  height = 6,
  width = 8,
)
