# Description: Creates folder structure for project.

library(purrr); library(fs); library(here)
purrr::walk(.x = fs::dir_ls('R'), .f = source)

# create directories for data and data-raw
dir_create(here("data", "genomic"))
dir_create(here("data-raw", "genomic"))
dir_create('data', 'dmet')
dir_create('data', 'drug_mapping')
dir_create('data', 'dmet', 'lines_of_therapy')


dir_create('output', 'report')
dir_create('output', 'fig')
dir_create('output', 'other')

fs::dir_create(here('data', 'survival'))
fs::dir_create(here('data', 'survival', 'first_line_platinum'))
fs::dir_create(here('data', 'survival', 'hrd_onco'))


