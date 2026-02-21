library(purrr)
library(fs)
library(here)
purrr::walk(.x = fs::dir_ls('R'), .f = source)

ca_ind <- readr::read_rds(
  here('data', 'cohort', 'ca_ind.rds')
)

met_time_custom <- readr::read_rds(
  file = here('data', 'cohort', 'met_time_custom.rds')
)

met_sites_at_met <- get_dmet_time(ca_ind_dat = ca_ind, annotate_type = T)
met_sites_at_met <- ca_ind %>%
  select(record_id, ca_seq, dob_ca_dx_days) %>%
  left_join(met_sites_at_met, ., by = c('record_id', 'ca_seq'))

met_sites_at_met %<>%
  mutate(
    dx_dmet_days = round(dx_dmet_yrs * 365.25),
    dob_dmet_days = dx_dmet_days + dob_ca_dx_days
  )


if (any(pull(count(met_sites_at_met, record_id)) > 1, na.rm = T)) {
  cli_abort('duplicate record_id in met_sites_at_met - fix!')
} else {
  met_sites_at_met %<>% select(-ca_seq)
}

met_sites_at_met <- met_time_custom %>%
  pivot_longer(
    cols = -record_id,
    names_to = 'site',
    values_to = 'dob_dmet_site_days'
  ) %>%
  left_join(
    met_sites_at_met,
    by = 'record_id'
  )

gg_dx_dmet_problem <- met_sites_at_met %>%
  mutate(diff = dob_dmet_site_days - dob_dmet_days) %>%
  filter(!is.na(diff)) %>%
  mutate(diff = diff + 1) %>%
  ggplot(aes(x = diff, y = site, color = .met_type)) +
  annotate(
    "rect",
    xmin = 0.5,
    xmax = 1.5,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.2,
    fill = "grey50"
  ) +
  geom_jitter(height = 0.3, width = 0.1, size = 0.5, alpha = 1) +
  scale_x_log10() +
  labs(
    x = "Days (log scale)",
    y = "Site (custom)",
    title = "Distribution of time from dmet (ca dx based) to site-specific met (img based)"
  ) +
  theme_minimal() +
  facet_wrap(vars(.met_type))

manu_plot_save_helper(
  dir = here('output', 'fig'),
  gg_dx_dmet_problem,
  name = 'gg_dx_dmet_problem',
  as_jpeg = F,
  as_png = F,
  height = 6,
  width = 8,
)

# We've shown why this doesn't really work.
# Our plan now is to use an imaging-based definition of first met.
met_sites_at_met %<>%
  # clear out the other defs:
  select(-c(dx_dmet_yrs, dob_ca_dx_days, dx_dmet_days, dob_dmet_days))


met_sites_at_met %<>%
  group_by(record_id) %>%
  mutate(
    dob_first_img_met_days = min(dob_dmet_site_days, na.rm = T)
  ) %>%
  ungroup(.)

dx_buffer = 0.5
binary_met_near_dx <- met_sites_at_met %>%
  mutate(
    has_met = case_when(
      is.na(dob_dmet_site_days) ~ F,
      dob_dmet_site_days <= (dob_first_img_met_days + dx_buffer) ~ T,
      T ~ F
    )
  ) %>%
  select(record_id, site, has_met, .met_type) %>%
  pivot_wider(
    names_from = site,
    values_from = has_met
  )

# binary_met_near_dx %>%
#   select(-c(record_id, .met_type))  %>%
#   gtsummary::tbl_summary(.)

binary_met_ever <- met_sites_at_met %>%
  mutate(
    has_met = case_when(
      is.na(dob_dmet_site_days) ~ F,
      T ~ T
    )
  ) %>%
  select(record_id, site, has_met) %>%
  pivot_wider(
    names_from = site,
    values_from = has_met
  )

binary_met_dual_col <- bind_rows(
  (binary_met_near_dx %>%
    mutate(epoch = '1st img date')),
  (binary_met_ever %>%
    mutate(epoch = 'Ever'))
) %>%
  mutate(epoch = forcats::fct_inorder(epoch))

binary_met_dual_col %>%
  select(-c(record_id, .met_type)) %>%
  gtsummary::tbl_summary(., by = epoch)
