# Do all the data manipulation needed to get a nice demo table.

library(purrr); library(fs); library(here);
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

read_wrap <- function(p) {
  read_rds(
    file = here("data", 'cohort', paste0(p, '.rds'))
  )
}

dft_pt <- read_wrap("pt")
dft_ca_ind <- read_wrap("ca_ind")
dft_med_onc <- read_wrap("med_onc")


dft_pt_baseline_sub <- dft_pt %>%
  mutate(
    `Race (primary)` = format_ptlevel_naaccr_race_code_primary(
      naaccr_race_code_primary
    ),
    `Ethnicity` = format_ptlevel_naaccr_ethnicity_code(
      naaccr_ethnicity_code,
    ),
    `Sex at birth` = factor(naaccr_sex_code)
  ) %>%
  select(record_id, 
         Institution = institution,
         `Race (primary)`,
         `Ethnicity`,
         `Sex at birth`,
         birth_year)

dft_ca_ind_baseline_sub <- dft_ca_ind %>%
  mutate(
    stage_mets_dx = format_stage_mets_dx(
      var_stage_dx = stage_dx,
      var_ca_dmets_yn = ca_dmets_yn
    ),
    ca_dx_how = format_ca_dx_how(ca_dx_how),
    ca_dmets_yn = format_ca_dmets_yn(ca_dmets_yn)
  ) %>%
  select(
    record_id,
    `Age at dx (years)` = dob_ca_dx_yrs,
    `Cancer type` = ca_type,
    `Cancer site (detailed)` = ca_d_site_txt,
    `Stage at dx` = stage_mets_dx
  ) 

# just forming a good spine for this:
dft_first <- dft_ca_ind %>%
  select(record_id) %>%
  distinct(.)

dft_first_ecog <- dft_med_onc %>%
  augment_med_onc_imputed_ecog(.) %>%
  filter(md_ecog_imputed != "Not documented in note") %>%
  group_by(record_id) %>%
  mutate(abs_days = abs(dx_md_visit_days)) %>%
  arrange(abs_days) %>% 
  slice(1) %>%
  ungroup(.) %>%
  select(
    record_id, 
    `First observed ECOG` =  md_ecog_imputed,
    # I'm leaving this in for now because I think it's interesting:
    `First ECOG source` = md_ecog_imp_source,
  )

dft_first %<>% left_join(., dft_first_ecog, by = 'record_id')

dft_met_ever <- get_dmet_time(dft_ca_ind) %>%
  group_by(record_id) %>%
  summarize(
    `Dmet anytime` = T
  )

dft_first %<>% left_join(., dft_met_ever, by = 'record_id')


# To find if MIBC was ever diagnosed we will use the "corrected" version,
# which considers either a med onc note or a stage >= 2 at diagnosis to be a sign on MIBC.
dft_mibc_ever <- make_dmet_musc_prop_status_block(dft_ca_ind) %>%
  group_by(record_id) %>%
  summarize(
    med_onc_mibc = any(status %in% c("Invasive", "Metastatic")),
    .groups = "drop"
  ) %>%
  left_join(
    ., select(dft_ca_ind_baseline_sub, record_id, `Stage at dx`), by = "record_id"
  ) %>%
  mutate(
    mibc_ever = case_when(
      med_onc_mibc ~ T,
      `Stage at dx` %in% c(
        "Stage II/III",
        "Stage IV (no met)",
        "Stage IV (met at dx)",
        "Stage IV (met unk.)"
      ) ~ T,
      T ~ F
    )
  ) %>%
  select(record_id, `MIBC anytime` = mibc_ever)

dft_first %<>% left_join(., dft_mibc_ever, by = 'record_id')

dft_first %<>%
  mutate(
    `First observed ECOG` = format_md_ecog(
      `First observed ECOG`,
      cluster_high_levels = T,
      na_lab = "(none)"
    ),
    `MIBC anytime` = case_when(
      is.na(`MIBC anytime`) ~ "No",
      `MIBC anytime` ~ "Yes", 
      T ~ "No"
    ),
    `Dmet anytime` = case_when(
      is.na(`Dmet anytime`) ~ "No",
      `Dmet anytime` ~ "Yes", 
      T ~ "No"
    )
  ) 






dft_demo <- full_join(
  dft_pt_baseline_sub,
  dft_ca_ind_baseline_sub,
  by = "record_id"
) %>%
  full_join(
    .,
    dft_first,
    by = "record_id"
  )



# age_dx is not an integer in this cohort, so this should be more exact.
dft_demo %<>% 
  mutate(
    `Year of birth` = birth_year,
    `Year of diagnosis` = round(birth_year + `Age at dx (years)`)
  ) %>%
  select(-birth_year) %>%
  relocate(`Year of birth`, `Year of diagnosis`,
           .before = `Age at dx (years)`) 


dft_demo %<>%
  relocate(
    `Age at dx (years)`, `Sex at birth`, 
    .before = 1L
  )

readr::write_rds(
  dft_demo,
  here('data', 'cohort', 'formatted_characteristics.rds')
)



# While working on changes (comment out otherwise):
dft_demo %>%
  select(-record_id) %>%
  gtsummary::tbl_summary(
    data = .,
    digits = list(
      `Year of birth` ~ 0,
      `Year of diagnosis` ~ 0
    )
  ) 

    

