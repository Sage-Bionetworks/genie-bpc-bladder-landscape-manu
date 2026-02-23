# We want to know how many NGS samples can be traced to a single anatomical location.
# This is somewhat complex actually in PRISSMM.  Each site can have multiple samples, but multiple samples could map to the same anatomic region.

library(purrr)
library(fs)
library(here)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

cpt <- readr::read_rds(
  here('data', 'cohort', 'cpt_aug.rds')
)

path <- readr::read_rds(
  here('data', 'cohort', 'path.rds')
)

path_sites_long <- path %>%
  select(
    cohort,
    record_id,
    path_proc_number,
    path_rep_number,
    contains('path_site')
  ) %>%
  pivot_longer(
    cols = contains('path_site'),
    names_to = 'path_site_num',
    values_to = 'icd_str'
  )

path_sites_long %<>%
  mutate(path_site_num = str_replace(path_site_num, "path_site", "")) %>%
  filter(!is.na(icd_str))

path_sites_long %<>%
  separate_wider_delim(
    cols = icd_str,
    delim = " ",
    names = c('icd_code', 'icd_desc'),
    too_many = 'merge',
    cols_remove = FALSE
  ) %>%
  mutate(
    icd_code = str_trim(icd_code),
    icd_desc = str_trim(icd_desc)
  )

# merge in our custom site categories:
icd_custom <- readr::read_csv(
  here('data-raw', 'manual', 'icd_map_custom.csv')
)

path_sites_long <- left_join(
  path_sites_long,
  select(icd_custom, icd_code, is_local, map_custom),
  by = 'icd_code',
  relationship = 'many-to-one'
)

# Ughhhhh.  There are 3 cases where a path proc/report has two associated ID numbers.
# For this reason we have to allow many-to-many here.
# It is also expected that the resulting data is much smaller than path_sites_long because tons of the path data has nothing to do with cpt samples, which sounds right.

cpt_path_long <-
  left_join(
    (cpt %>%
      select(
        cpt_genie_sample_id,
        cohort,
        record_id,
        path_proc_number,
        path_rep_number
      )),
    path_sites_long,
    by = c('cohort', 'record_id', 'path_proc_number', 'path_rep_number'),
    relationship = 'many-to-many'
  )

readr::write_rds(
  cpt_path_long,
  here('data', 'cohort', 'cpt_path_long.rds')
)
