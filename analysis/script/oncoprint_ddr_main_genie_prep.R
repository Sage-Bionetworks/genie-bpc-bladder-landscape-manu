library(purrr)
library(here)
library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)
# little weird - we're pulling the first samples because that file has the record_id linked in.
blad_samp_mg_first <- readr::read_rds(
  here('data', 'genomic', 'main_genie', 'bladder_samples_mg_first.rds')
)

writeLines(
  text = unique(blad_samp_mg_first$record_id),
  con = here('data', 'genomic', 'oncoprint', 'ddr_mg', 'oncoprint_pt.txt')
)

# Just going to pull the gene file from the BPC DDR oncoprint folder
fs::file_copy(
  path = here(
    'data',
    'genomic',
    'oncoprint',
    'ddr',
    'oncoprint_genes.txt'
  ),
  new_path = here(
    'data',
    'genomic',
    'oncoprint',
    'ddr_mg',
    'oncoprint_genes.txt'
  )
)
