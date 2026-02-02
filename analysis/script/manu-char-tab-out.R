library(purrr)
library(fs)
library(here)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

theme_gtsummary_language("en", big.mark = "")

demo_coh <- readr::read_rds(
  here('data', 'cohort', 'formatted_characteristics.rds')
)

demo_coh %>%
  select(
    -c(
      record_id,
      # taking off some stuff that's too detailed:
      `Cancer site (detailed)`,
      `First ECOG source`
    )
  ) %>%
  gtsummary::tbl_summary(
    data = .,
    digits = list(
      `Year of birth` ~ 0,
      `Year of diagnosis` ~ 0
    )
  ) %>%
  as_gt(.) %>%
  gt::gtsave(
    filename = here('output', 'manu', 'supp', 'SX_char_table_cohort.docx')
  )

demo_mg <- readr::read_rds(
  here('data', 'cohort', 'formatted_characteristics_with_main_genie.rds')
)

demo_mg %>%
  select(study, `Age at seq.`:Institution) %>%
  gtsummary::tbl_summary(
    data = .,
    by = study
  ) %>%
  as_gt %>%
  gt::gtsave(
    filename = here(
      'output',
      'manu',
      'supp',
      'SX_char_table_vs_main_genie.docx'
    )
  )
