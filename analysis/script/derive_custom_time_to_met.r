
library(purrr); library(fs); library(here)
purrr::walk(.x = fs::dir_ls('R'), .f = source)

icd_custom <- readr::read_csv(
  here('data-raw', 'manual', 'icd_map_custom.csv')
)

img <- readr::read_rds(
  here('data', 'cohort', 'img.rds')
)
ca_ind <- readr::read_rds(
  here('data', 'cohort', 'ca_ind.rds')
)

img_ca <- img %>%
  # not using dx_scan_days because its from the first bpc project cancer.
  #   can't imagine a scenario where that's useful.
  select(record_id, scan_number, image_scan_int, image_ca, 
         matches('image_casite[0-9]{1,2}'))

img_ca %<>%
  pivot_longer(
    cols = matches('image_casite[0-9]{1,2}'),
    names_to = 'loc_num',
    values_to = 'icd_str'
  ) %>%
  filter(!is.na(icd_str)) 

# varying amount of whitespace in here, so we'll split and trim in steps.
img_ca %<>% 
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

img_ca <- left_join(
  img_ca,
  select(icd_custom, icd_code, is_local, map_custom),
  by = "icd_code"
)


# for each person find the time to each met category.
# we can always make a category for local lymph if needed, for but now those
#  will be removed.
img_ca %<>% filter(!is_local)

met_time_custom <- img_ca %>% 
  group_by(record_id, map_custom) %>%
  summarize(
    dob_met_days = min(image_scan_int, na.rm = T),
    .groups = 'drop'
  )

# across everything, find the people with mets.
met_time_total <- met_time_custom %>%
  group_by(record_id) %>%
  summarize(
    dob_met_days = min(dob_met_days, na.rm = T)
  )

diff_derived_custom <- setdiff(
  get_dmet_time(ca_ind)$record_id, met_time_total$record_id
)
diff_custom_derived <- setdiff(
  met_time_total$record_id, get_dmet_time(ca_ind)$record_id
)

ca_ind %>%
  filter(record_id %in% diff_derived_custom) %>%
  count(stage_dx_iv)

ca_ind %>%
  filter(record_id %in% diff_derived_custom) %>%
  select(record_id, matches('dx_to_dmets_.*_days')) %>%
  mutate(
    across(
      .cols = matches('dx_to_dmets_.*_days'),
      .fns = \(x) !is.na(x)
    )
  ) %>%
  pivot_longer(
    cols = -record_id
  ) %>%
  filter(value) %>% # only those with a met in that location
  count(name) # number of people with mets at each spot
# abdomen we know about - renal pelvis is the big one I know of.
# Pelvis is pretty vague and not interesting, so we can probably move on there.


met_time_custom <- met_time_custom %>%
  mutate(
    map_custom = case_when(
      map_custom %in% "Lymph node (distant)" ~ "lymph_distant",
      T ~ tolower(map_custom)
    )
  ) %>%
  pivot_wider(
    names_from = 'map_custom',
    values_from = 'dob_met_days'
  )

readr::write_rds(
  x = met_time_custom,
  file = here('data', 'cohort', 'met_time_custom.rds')
)
