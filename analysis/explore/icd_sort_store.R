library(here)
library(readr)
library(dplyr)
library(tidyr)

icd_map <- readr::read_csv(
  here('analysis', 'explore', 'bpc_bladder_icd_map.csv')
)

icd_map <- icd_map |>
  replace_na(list(is_local = FALSE)) |>
  arrange(is_local, map_custom, icd_code) |>
  select(is_local, map_custom, everything())

readr::write_csv(
  icd_map,
  here('data-raw', 'manual', 'icd_map_custom.csv'),
  na = ""
)
