
library(purrr); library(fs); library(here)
purrr::walk(.x = fs::dir_ls('R'), .f = source)

map_custom <- readxl::read_xls(
  here('data-raw', 'manual', 'New_Mappings_njs.xls')
)

map_custom %<>% 
  select(Kode, Map) %>%
  distinct(.)

map_custom %>%
  filter(Kode %in% map_custom$Kode[duplicated(map_custom$Kode)]) %>%
  arrange(Kode)

# temporary fix for the above
map_custom %<>%
  arrange(Kode, Map %in% "No map") %>%
  group_by(Kode) %>%
  slice(1) %>%
  ungroup(.)

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

# Just a check to make sure I understand the structure of this data:
img_ca %>% 
  filter(image_ca %in% "Yes, the Impression states or implies there is evidence of cancer")

# These are not expected to me - there is cancer but no sites are mentioned as having cancer.  I guess cancer noted with no known ICD code?
img_ca %>% 
  filter(
    image_ca %in% "Yes, the Impression states or implies there is evidence of cancer" &
      is.na(image_casite1)
  ) %>% glimpse

# All the same, I guess we can proceed for the purposes of identifying met sites.

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

img_ca %<>% 
  left_join(
    .,
    map_custom,
    by = c(icd_code = "Kode")
  )

parsimonious_mapping_file <- img_ca %>% 
  count(icd_code, icd_desc, name = 'n_times_reported') %>%
  left_join(
    map_custom,
    by = c(icd_code = "Kode")
  ) %>%
  arrange(Map)

readr::write_csv(
  parsimonious_mapping_file,
  here('analysis', 'explore', 'bpc_bladder_icd_map.csv')
)

img_ca %>% count(Map, sort = T)

img_ca %>% 
  mutate(
    Map = case_when(
      str_detect(Map, 'Lymph node') ~ "lymph node",
      T ~ Map
    )
  ) %>%
  group_by(record_id, Map) %>%
  summarize(any_met = TRUE, .groups = 'drop') %>%
  count(Map) %>%
  print(n = 500)


img_ca %>% 
  mutate(
    Map = case_when(
      str_detect(Map, 'Lymph node') ~ "lymph node",
      T ~ Map
    )
  ) %>%
  group_by(record_id) %>%
  summarize(only_lymph = , .groups = 'drop') %>%
  count(Map) %>%
  print(n = 500)

ca_ind %>%
  summarize(
    across(
      .cols = c(dmets_brain, dmets_bone, dmets_liver, dmets_thorax),
      .fns = sum
    )
  )

