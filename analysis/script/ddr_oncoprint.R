# The top-genes oncoprint is imbedded within a script because it's
#   somewhat complicated to find the top 10 most altered genes.
# This one can be a standalone because our definiton of DDR is easy.

library(purrr)
library(here)
library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

read_wrap <- function(p) {
  read_csv(file = here("data-raw", p), show_col_types = F)
}

pt <- read_rds(
  file = here("data", 'cohort', 'pt.rds')
)

writeLines(
  text = unique(pt$record_id), # should be unique already
  con = here('data', 'genomic', 'oncoprint', 'ddr', 'oncoprint_pt.txt')
)

# Here I'm choosing the order to give % sorting on the oncoprint.
ddr_9_panel_genes <- c(
  'ATM',
  'BRCA2',
  'ATR',
  'BRCA1',
  'ERCC2',
  'ERCC5',
  'RAD51C',
  'FANCC',
  'RECQL4'
)

writeLines(
  text = ddr_9_panel_genes,
  con = here('data', 'genomic', 'oncoprint', 'ddr', 'oncoprint_genes.txt')
)
