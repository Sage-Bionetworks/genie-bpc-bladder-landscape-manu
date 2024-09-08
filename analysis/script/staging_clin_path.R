library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

dft_ca_ind <- readr::read_rds(here('data', 'cohort', "ca_ind.rds"))

# We need a TNM grouping at both the clinical and pathological staging steps
#  for the analysis we have planned.  First step is to check through the ca_ind
#  variables to see if we have this already coded somehow.

# Candidate variables:
# - best_ajcc_stage_cd 
# - stage_dx:  nope, combines everything.
# - ca_tnm_edition_num:  theoretically useful, I'm just going use 8.
# - naaccr_tnm_path_desc:  No, this is stage not TNM categories.
# - naaccr_clin_stage_cd:  This looks great, except we only have the clinical stage.
# So no, we don't already ahve this.  Great, let's get to work.

dft_tnm <- dft_ca_ind %>% 
  select(
    record_id,
    ca_seq,
    naaccr_path_t_cd,
    ca_path_t_stage,
    naaccr_clin_t_cd,
    ca_clin_t_stage,
    
    naaccr_path_n_cd,
    ca_path_n_stage,
    naaccr_clin_n_cd,
    ca_clin_n_stage,
    
    # We don't seem to have the ca_*_m_stage variables for M categories,
    #   have to make due with the tumor registry stuff I guess.
    naaccr_clin_m_cd,
    naaccr_path_m_cd
  )

dft_tnm_before <- dft_tnm

dft_tnm %<>%
  mutate(
    across(
      .cols = -c(record_id, ca_seq),
      .fns = \(z) if_else(
        z %in% c("Unknown", "Not Applicable"),
        NA_character_,
        z
      )
    )
  ) %>%
  mutate(
    across(
      .cols = -c(record_id, ca_seq),
      .fns = \(z) case_when(
        str_detect(z, "X") ~  NA_character_,
        # for 88 here's my source:
        # https://staging.seer.cancer.gov/naaccr/item/tnm/1.9/940/
        str_detect(z, "88") ~  NA_character_,
        T ~ z
      )
    )
  )

# We're going to use a function to strip out all the stuff that's in front.
# This is because we have rediculous gibberish entries like "NcN0" and "N88"
#   which have been added by layers of programming without error checks, I assume.
# The stuff we'll strip off is C, clinical, P, pathological (maybe?), T,N and M
tnm_stripper <- function(vec) {
  # The stuff we'll strip off is C, clinical, P, pathological (maybe?), T,N and M.
  # yc and yp are used for post-therapy and post-neoadjuvant, which probably
  #   shoudn't be here but that's fine.  rc and rp are used for recurrence,
  #   again shouldn't be here.
  str_replace_all(
    vec, "^[ryPpCcTtNnMm]*", ""
  )
}


dft_tnm %<>%
  mutate(
    across(
      .cols = -c(record_id, ca_seq),
      .fns = tnm_stripper
    )
  ) %>%
  # at this point if we have a blank string, they wrote something like "N",
  #   which is less helpful than writing nothing.
  mutate(
    across(
      .cols = -c(record_id, ca_seq),
      .fns = \(z) case_when(
        z %in% "" ~  NA_character_,
        # also taking this chance to lowercase everything now.
        # There are no valid capitals AFTER the number, so this is just righting wrongs.
        T ~ tolower(z)
      )
    )
  )


dft_tnm_after_stripping <- dft_tnm

# QC step: look at how the codes changed under these transformations:
dft_tnm_qc <- full_join(
  pivot_longer(
    dft_tnm_before,
    cols = -c(record_id, ca_seq),
    values_to = "before"
  ),
  pivot_longer(
    dft_tnm_after_stripping,
    cols = -c(record_id, ca_seq),
    values_to = "after"
  ),
  by = c('record_id', 'ca_seq', 'name')
) %>%
  count(name, before, after, name = "cases") %>%
  select(cases, everything()) %>%
  arrange(name, after, before)

readr::write_rds(
  dft_tnm_qc,
  here('data', 'cohort', 'tnm_path_clin_stripping_qc.rds')
)

# prioritizing the tumor registry over curated variables just follows the way
#   the curation is done - they skip it if the tumor registry is filled out.
dft_tnm %<>%
  mutate(
    clin_comb_n = case_when(
      !is.na(naaccr_clin_n_cd) ~ naaccr_clin_n_cd,
      T ~ ca_clin_n_stage
    ),
    path_comb_n = case_when(
      !is.na(naaccr_path_n_cd) ~ naaccr_path_n_cd,
      T ~ ca_path_n_stage
    ),
    clin_comb_t = case_when(
      !is.na(naaccr_clin_t_cd) ~ naaccr_clin_t_cd,
      T ~ ca_clin_t_stage
    ),
    path_comb_t = case_when(
      !is.na(naaccr_path_t_cd) ~ naaccr_path_t_cd,
      T ~ ca_path_t_stage
    ),
    # only one game in town here:
    clin_comb_m = naaccr_clin_m_cd,
    path_comb_m = naaccr_path_m_cd
  )

# these sets of valid values come from the AJCC manual, 8th edition, starting on
#   page 771 ("definition of primary tumor (T)")
ajcc_bladder_t_cat <- c("TX", "T0", "Ta", "Tis", "T1", "T2",
                        "pT2a", "pT2b", "T3", "pT3a", "pT3b",
                        "T4", "T4a", "T4b")
s_ajcc_bladder_t_cat <- tnm_stripper(ajcc_bladder_t_cat)
ajcc_bladder_n_cat <- c("NX", "N0", "N1", "N2", "N3")
ajcc_bladder_m_cat <- c("M0", "M1", "M1a", "M1b")
# Probably need to go back and check how the renal pelvis cases line up.

lev_stage_groups <- c(
  # I added zero here, because we have to do something with people who are 
  # T0N0M0 (quite a few).
  "0", "0a", "0is", "I", "II", 'IIIA', 'IIIB', 'IVA', 'IVB'
)

dft_tnm %<>%
  mutate(
    # This comes from the AJCC on page 772, "AJCC prognostic stage groups"
    clin_group_stage = case_when(
      clin_comb_m %in% "1b" ~ "IVB",
      clin_comb_m %in% "1a" ~ "IVA",
      clin_comb_m %in% "1" ~ "IVA", # unclear, going with the lower.
      clin_comb_t %in% "4b" ~ "IVA",
      # for all of the remaining categories I'm making the choice to not check M.
      # this is due to the fact that silence is a zero for M in the latest editions,
      #  which differs from T and N.
      clin_comb_t %in% s_ajcc_bladder_t_cat[5:13] &
        clin_comb_n %in% c("2", "3") ~ "IIIB",
      clin_comb_t %in% s_ajcc_bladder_t_cat[5:13] &
        clin_comb_n %in% c("1") ~ "IIIA",
      clin_comb_t %in% s_ajcc_bladder_t_cat[9:13] &
        clin_comb_n %in% c("0") ~ "IIIA",
      clin_comb_t %in% s_ajcc_bladder_t_cat[6:8] &
        clin_comb_n %in% c("0") ~ "II",
      clin_comb_t %in% "1" & clin_comb_n %in% c("0") ~ "I",
      clin_comb_t %in% "is" & clin_comb_n %in% c("0") ~ "0is",
      clin_comb_t %in% "a" & clin_comb_n %in% c("0") ~ "0a",
      clin_comb_t %in% "0" & clin_comb_n %in% c("0") ~ "0",
      T ~ NA_character_
    )
  )

# 100% copypasta from above.  
dft_tnm %<>%
  mutate(
    # This comes from the AJCC on page 772, "AJCC prognostic stage groups"
    path_group_stage = case_when(
      path_comb_m %in% "1b" ~ "IVB",
      path_comb_m %in% "1a" ~ "IVA",
      path_comb_m %in% "1" ~ "IVA", # unclear, going with the lower.
      path_comb_t %in% "4b" ~ "IVA",
      # for all of the remaining categories I'm making the choice to not check M.
      # this is due to the fact that silence is a zero for M in the latest editions,
      #  which differs from T and N.
      path_comb_t %in% s_ajcc_bladder_t_cat[5:13] &
        path_comb_n %in% c("2", "3") ~ "IIIB",
      path_comb_t %in% s_ajcc_bladder_t_cat[5:13] &
        path_comb_n %in% c("1") ~ "IIIA",
      path_comb_t %in% s_ajcc_bladder_t_cat[9:13] &
        path_comb_n %in% c("0") ~ "IIIA",
      path_comb_t %in% s_ajcc_bladder_t_cat[6:8] &
        path_comb_n %in% c("0") ~ "II",
      path_comb_t %in% "1" & path_comb_n %in% c("0") ~ "I",
      path_comb_t %in% "is" & path_comb_n %in% c("0") ~ "0is",
      path_comb_t %in% "a" & path_comb_n %in% c("0") ~ "0a",
      path_comb_t %in% "0" & path_comb_n %in% c("0") ~ "0",
      T ~ NA_character_
    )
  )

lev_stage_clust <- c("0/I", "II", "IIIA", "IIIB", "IVA", "IVB")

dft_tnm %<>%
  mutate(
    clin_group_clust = case_when(
      is.na(clin_group_stage) ~ NA_character_,
      clin_group_stage %in% lev_stage_groups[1:4] ~ "0/I",
      T ~ clin_group_stage
    ),
    path_group_clust = case_when(
      is.na(path_group_stage) ~ NA_character_,
      path_group_stage %in% lev_stage_groups[1:3] ~ "0/I",
      T ~ path_group_stage
    )
  )
      

dft_tnm %<>%
  mutate(
    across(
      .cols = c(clin_group_stage, path_group_stage),
      .fns = \(z) factor(z, levels = lev_stage_groups, ordered = T)
    ),
    across(
      .cols = c(clin_group_clust, path_group_clust),
      .fns = \(z) factor(z, levels = lev_stage_clust, ordered = T)
    )
  )


dft_tnm_classification_qc <- dft_tnm %>% 
  select(
    record_id, ca_seq,
    contains('_comb_'),
    contains('group_clust')
  ) %>%
  mutate(across(-ca_seq, as.character)) %>%
  pivot_longer(
    cols = -c(record_id, ca_seq)
  ) %>%
  separate_wider_delim(
    cols = name,
    delim = "_",
    names = c("time", "obs"),
    too_many = "merge"
  ) %>%
  pivot_wider(
    names_from = 'obs',
    values_from = 'value'
  ) %>%
  count(time, comb_t, comb_n, comb_m, group_clust, name = "cases") %>%
  select(cases, everything()) %>%
  arrange(group_clust, comb_m, comb_n, comb_t)

dft_tnm_classification_qc %>% print(n = 500)
    

readr::write_rds(
  dft_tnm,
  here('data', 'cohort', 'tnm_path_clin_details.rds')
)
  



