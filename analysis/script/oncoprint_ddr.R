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

ddr_9_panel_genes <- c(
  'ATM',
  'ATR',
  'BRCA1',
  'BRCA2',
  'RECQL4',
  'RAD51C',
  'FANCC',
  'ERCC2',
  'ERCC5'
)
