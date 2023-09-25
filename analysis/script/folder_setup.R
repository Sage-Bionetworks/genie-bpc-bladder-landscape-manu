# Description: Creates folder structure for project.

library(purrr); library(fs); library(here)
purrr::walk(.x = fs::dir_ls('R'), .f = source)

# create directories for data and data-raw
dir_create(here("data", "genomic"))
dir_create(here("data-raw", "genomic"))
dir_create('data', 'dmet')

