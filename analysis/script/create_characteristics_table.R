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
# dft_img <- read_wrap("img")
dft_med_onc <- read_wrap("med_onc")
# dft_path <- read_wrap("path")
# dft_reg <- read_wrap("reg")
# dft_tm <- read_wrap("tm_level_dataset.csv")
# dft_cpt <- read_wrap("cpt")
# dft_drug <- read_wrap("drug")
# dft_rad <- read_wrap("rad")


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
    mets_s4 = s4_met_helper(stage_dx_iv, ca_dmets_yn),
    ca_dx_how = format_ca_dx_how(ca_dx_how),
    ca_hist_adeno_squamous = format_ca_hist_adeno_squamous(ca_hist_adeno_squamous),
    stage_dx = format_stage_dx(stage_dx),
    ca_dmets_yn = format_ca_dmets_yn(ca_dmets_yn),
    ca_path_n_stage = format_ca_path_n_stage(ca_path_n_stage),
  ) %>%
  select(
    record_id,
    `Age at dx (years)` = dob_ca_dx_yrs,
    `Cancer type` = ca_type,
    `Stage at dx` = stage_dx,
    # Removing this:
    # `Source of dx date` = ca_dx_how,
    `Mets at dx (for Stage IV pts.)` = mets_s4
    # Update May 5:  Removed histology (unimportant)
    # `Histology` = ca_hist_adeno_squamous, # at dx?
    # Update Mar 23 2024:  remove this too (unimportant)
    # `Pathologic N Stage` = ca_path_n_stage # describes spread to lymph nodes
  ) 


dft_med_onc_dx <- augment_med_onc_imputed_ecog(dft_med_onc)

dft_med_onc_dx %<>%
  filter(md_ecog_imputed != "Not documented in note") %>%
  mutate(md_ecog_imputed = format_md_ecog(md_ecog_imputed)) %>%
  group_by(record_id) %>%
  filter(dx_md_visit_days < 30 & dx_md_visit_days > -180) %>%
  mutate(abs_days = abs(dx_md_visit_days)) %>%
  arrange(abs_days) %>% 
  slice(1) %>%
  ungroup(.) %>%
  select(
    record_id, 
    `ECOG (or Karnofsky)` =  md_ecog_imputed,
    `ECOG scale source` = md_ecog_imp_source,
  )

dft_demo <- full_join(
  dft_pt_baseline_sub,
  dft_ca_ind_baseline_sub,
  by = "record_id"
) %>%
  full_join(
    .,
    dft_med_onc_dx,
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

readr::write_rds(
  dft_demo
  here('data', 'cohort', 'formatted_characteristics.rds')
)