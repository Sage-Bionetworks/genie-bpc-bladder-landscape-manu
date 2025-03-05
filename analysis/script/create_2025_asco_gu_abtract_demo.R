library(purrr); library(fs); library(here);
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

ddr_flags <- readr::read_rds(
  here('data', 'genomic', 'ddr_def_compare', 'ddr_flags_first_sample.rds')
)

ddr_flags %<>% rename(record_id = patient_id)

cohort_char <- readr::read_rds(
  here('data', 'cohort', 'formatted_characteristics.rds')
)
cohort_char %<>% select(
  -c(`Cancer site (detailed)`,
     `First ECOG source`)
)

met_time_custom <- readr::read_rds(
  file = here('data', 'cohort', 'met_time_custom.rds')
)
met_flag <- met_time_custom %>%
  # just the ones we happen to have on the poster
  select(record_id, lymph_distant, 
         bone, lung, liver) %>%
  mutate(
    across(.cols = -record_id,
           .fns = \(x) if_else(is.na(x), F, T)
    )
  ) %>%
  rename_with(~str_c("Met_anytime_", .), -record_id)

ddr_demo <- ddr_flags %>%
  left_join(cohort_char, by = 'record_id') %>%
  left_join(met_flag, by = 'record_id')

ddr_demo %<>%
  mutate(
    across(
      matches("^Met_"),
      .fns = \(x) if_else(is.na(x), F, x)
    )
  )

ddr_demo <- bind_rows(
  (ddr_demo %>%
     mutate(ddr = if_else(ddr, "DDR+", "DDR-"))),
  (ddr_demo %>% mutate(ddr = "All"))
) %>%
  mutate(ddr = fct_inorder(ddr))


theme_gtsummary_compact(font_size = 12)
theme_gtsummary_language("en", big.mark = "") 

disp_demo_by_ddr <- gtsummary::tbl_summary(
  select(ddr_demo, -c('record_id', 'sample_id')),
  by = ddr,
  digits = list(
    `Year of birth` ~ 0,
    `Year of diagnosis` ~ 0
  )
) 


style_help <- function(ft) {
  ft %>% 
    theme_box(.) %>%
    fontsize(size = 6, part = "all") %>%
    valign(part = "all", valign = "top") %>%
    padding(part = 'body', padding = 1) %>% 
    merge_v(j = c(1)) %>%
    autofit(.) 
}

disp_demo_by_ddr %>%
  as_flex_table(.) %>%
  flextable::fontsize(size = 8, part = "all") %>%
  list(`Characteristics by DDR status of first sample` = .) %>%
  save_as_docx(
    values = .,
    path = here('output', 'other', 'asco_gu_2025_ddr_demo.docx')
  )
