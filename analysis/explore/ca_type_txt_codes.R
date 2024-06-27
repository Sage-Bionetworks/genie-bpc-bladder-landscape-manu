library(purrr); library(fs); library(here);
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

read_wrap <- function(p) {
  read_rds(
    file = here("data", 'cohort', paste0(p, '.rds'))
  )
}

dft_ca_ind <- read_wrap("ca_ind")

openxlsx::write.xlsx(
  x = count(dft_ca_ind, ca_d_site_txt, sort = T),
  file = here('analysis', 'explore', 'ca_d_site_codes.xlsx')
)

