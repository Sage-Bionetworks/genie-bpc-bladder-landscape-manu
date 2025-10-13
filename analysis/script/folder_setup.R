# Description: Creates folder structure for project.

library(purrr)
library(fs)
library(here)
purrr::walk(.x = fs::dir_ls('R'), .f = source)

# create directories for data and data-raw
dir_create('data', 'genomic', 'main_genie')
dir_create('data-raw', 'genomic', 'main_genie')
dir_create('data', 'dmet')
dir_create('data', 'drug_mapping')
dir_create('data', 'dmet', 'lines_of_therapy')
dir_create('data', 'genomic', 'gene_corr')

dir_create('output', 'report')
dir_create('output', 'fig')
dir_create('output', 'other')
dir_create('output', 'aacr_ss24', 'tfl_obj')
dir_create('output', 'aacr_ss24', 'img')

fs::dir_create(here('data', 'survival'))
fs::dir_create(here('data', 'survival', 'first_line_platinum'))
fs::dir_create(here('data', 'survival', 'hrd_onco'))
fs::dir_create(here('data', 'survival', 'ercc3_plat'))
fs::dir_create(here('data', 'survival', 'first_line_immuno'))
fs::dir_create(here('data', 'survival', 'ddr_onco'))
fs::dir_create(here('data', 'survival', 'ddr_neoadj'))
fs::dir_create(here('data', 'survival', 'line23'))

fs::dir_create(here('data', 'genomic', 'ddr_def_compare', 'ddr_as_outcome'))
