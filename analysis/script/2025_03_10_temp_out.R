library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

dir_out <- here('output', 'fig', '2025_03_10_temp')

bundle_1L <- readr::read_rds(
  here('data', 'survival', 'ddr_onco_1L', 'model_interact_bundle.rds')
)
bundle_2L <- readr::read_rds(
  here('data', 'survival', 'ddr_onco_2L', 'model_interact_bundle.rds')
)


gg_km_1L <- bundle_1L$gg_km_1L

ggsave(
  plot = bundle_1L$gg_km_1L,
  filename = here(dir_out, '1L_KM.pdf'),
  height = 6, width = 6
)
ggsave(
  plot = bundle_1L$gg_imputed,
  filename = here(dir_out, '1L_model.pdf'),
  height = 6, width = 6
)
ggsave(
  bundle_2L$gg_km_2L,
  filename = here(dir_out, '2L_KM.pdf'),
  height = 6, width = 6
)
ggsave(
  bundle_2L$gg_imputed,
  filename = here(dir_out, '2L_model.pdf'),
  height = 6, width = 6
)

ggsave(
  bundle_2L$gg_km_2L,
  filename = here(dir_out, '2L_KM.pdf'),
  height = 6, width = 6
)
ggsave(
  bundle_2L$gg_imputed,
  filename = here(dir_out, '2L_model.pdf'),
  height = 6, width = 6
)
