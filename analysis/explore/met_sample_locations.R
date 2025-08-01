cpt_aug <- readr::read_rds(
  here('data', 'cohort', 'cpt_aug.rds')
)

path <- readr::read_rds(
  here('data', 'cohort', 'path.rds')
)

# link between these two is record, path_proc_number, path_report_number

# bit of a problme here because if there are multiple anatomic locations that's
#   only linkable to a single cancer panel test - we don't know which one.
#   Perhaps they use a blender?

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
    values_to = 'icd_code'
  )

path_sites_long %<>%
  mutate(path_site_num = str_replace(path_site_num, "path_site", "")) %>%
  filter(!is.na(icd_code)) %>%
  # group by everything but path_site_num - drop multi-site things with the same anatomic location.
  group_by(cohort, record_id, path_proc_number, path_rep_number, icd_code) %>%
  arrange(path_site_num) %>%
  slice(1) %>%
  ungroup()
# actually drops a lot.

# Gives a quick and dirty "how many things with multiple sites do we have?"
nrow(path_sites_long)
path_sites_long %>%
  group_by(record_id, path_proc_number, path_rep_number) %>%
  count() %>%
  ungroup() %>%
  count(n) # about 50% are solo.
# Try again with urine excluded:
path_sites_long_no_urine <- path_sites_long %>%
  filter(!str_detect(icd_code, "F40"))
nrow(path_sites_long_no_urine)
path_sites_long_no_urine %>%
  group_by(record_id, path_proc_number, path_rep_number) %>%
  count() %>%
  ungroup() %>%
  count(n) # 35% are solo anatomic locations.

path_sites_long %>%
  count(icd_code, sort = T) %>%
  readr::write_csv(
    .,
    here('analysis', 'explore', 'path_report_locations.csv')
  )
