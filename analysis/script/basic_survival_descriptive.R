# Create some basic plots, KM curves, truncation tests, etc.
# Don't want these to clutter up the report too much.

library(purrr)
library(here)
library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source) # also loads lots of packages.

fs::dir_create(
  here('data', 'survival')
)

surv_desc_fp <- here('data', 'survival')

read_wrap_clin <- function(p) {
  read_rds(file = here("data", 'cohort', p))
}

dft_pt <- read_wrap_clin("pt.rds")
dft_ca_ind <- read_wrap_clin("ca_ind.rds")
# we're taking the augmented version, which has TMB columns added.  The existing info is the same.
dft_cpt <- read_wrap_clin("cpt_aug.rds")


# Create plot for survival from diagnosis stratified by stage at dx.
dft_first_cpt <- get_first_cpt(
  ca_ind_dat = dft_ca_ind,
  cpt_dat = dft_cpt
)

dft_surv_dx <- dft_ca_ind %>%
  left_join(., dft_first_cpt, by = c("record_id", "ca_seq")) %>%
  mutate(
    stage_custom = case_when(
      stage_dx %in% c("Stage 0", "Stage I") ~ "Early (Stage 0-1)",
      stage_dx %in% c("Stage II", "Stage III") ~ "Invasive (Stage 2-3)",
      stage_dx %in% c("Stage IV") ~ "Stage IV",
      T ~ NA_character_ # Note:  Stage I-III NOS and the 3 NAs get clumped here.
    ),
    stage_custom = factor(stage_custom) # Happens to be alphabetical.
  ) %>%
  select(
    record_id,
    ca_seq,
    stage_custom,
    dx_cpt_rep_yrs,
    os_dx_status,
    tt_os_dx_yrs
  )

dft_surv_dx %<>%
  remove_trunc_gte_event(
    trunc_var = 'dx_cpt_rep_yrs',
    event_var = 'tt_os_dx_yrs'
  )

surv_obj_os_dx <- with(
  dft_surv_dx,
  Surv(
    time = dx_cpt_rep_yrs,
    time2 = tt_os_dx_yrs,
    event = os_dx_status
  )
)

gg_os_dx_stage <- plot_one_survfit(
  dat = dft_surv_dx,
  surv_form = surv_obj_os_dx ~ stage_custom,
  plot_title = "OS from diagnosis",
  plot_subtitle = "Adjusted for (independent) delayed entry",
  x_exp = 0.1
)

readr::write_rds(
  x = gg_os_dx_stage,
  file = here(surv_desc_fp, 'os_dx_by_stage.rds')
)

gg_os_dx_stage_manu <- plot_one_survfit(
  dat = dft_surv_dx,
  surv_form = surv_obj_os_dx ~ stage_custom,
  plot_title = "OS from diagnosis",
  plot_subtitle = "Adjusted for (independent) delayed entry",
  x_exp = 0.1,
  x_breaks = seq(0, 5, by = 0.5),
  risktable_prop = 0.3,
  risktable_font_size = 3
) +
  coord_cartesian(xlim = c(0, 5), ylim = c(0, 1.01))

gg_os_dx_stage_manu <- gg_os_dx_stage_manu +
  labs(subtitle = NULL, title = "Overall survival from dx") +
  theme(plot.title.position = 'panel')

readr::write_rds(
  x = gg_os_dx_stage_manu,
  file = here(surv_desc_fp, 'os_dx_by_stage_manu.rds')
)

ggsave(
  gg_os_dx_stage_manu, # changed the name to "manu" as the project evolved.
  height = 5,
  width = 8,
  filename = here(
    'output',
    'aacr_ss24',
    'img',
    '03_surv_dx_stage.pdf'
  )
)


# Create survival variables from met.
# In prostate cancer this is different from advanced disease.

dft_surv_dmet <- left_join(
  dft_ca_ind,
  get_dmet_time(dft_ca_ind, annotate_type = T),
  by = c("record_id", "ca_seq")
) %>%
  filter(!is.na(dx_dmet_yrs)) # makes no sense otherwise.

dft_surv_dmet %<>%
  mutate(
    tt_os_dmet_days = tt_os_dx_days - (dx_dmet_yrs * 365.25),
    tt_os_dmet_mos = tt_os_dx_mos - (dx_dmet_yrs * 12.0148),
    tt_os_dmet_yrs = tt_os_dx_yrs - (dx_dmet_yrs)
  )

# Difference mean plot
ggplot(
  dft_surv_dmet,
  aes(
    y = tt_os_adv_days - tt_os_dmet_days,
    x = (tt_os_adv_days + tt_os_dmet_days) / 2,
    color = .met_type
  )
) +
  geom_jitter(alpha = 0.5, width = 0, height = 100) +
  theme_bw() +
  labs(
    title = "Diff/mean for 'advanced' vs 'metatstatic'",
    subtitle = "QA: Y ~= 0 for all but stage IV, dmet later"
  )

dft_surv_dmet %<>% select(-.met_type)

dft_surv_dmet %<>%
  left_join(., dft_first_cpt, by = c("record_id", "ca_seq")) %>%
  mutate(
    dmet_cpt_rep_yrs = dx_cpt_rep_yrs - dx_dmet_yrs
  )

## Heuristic:
# dft_surv_dmet %>%
#   arrange(desc(dmet_cpt_rep_yrs)) %>%
#   slice(1:3) %>%
#   select(record_id, contains("cpt"), dx_dmet_yrs, tt_os_dmet_yrs, os_dx_status) %>%
#   View(.)

readr::write_rds(
  x = dft_surv_dmet,
  file = here(surv_desc_fp, 'ca_ind_os_dmet.rds')
)


# Create a survival plot showing the effect of adjusting for truncation.
# We will go from metastasis.

pal_surv_dmet <- c('#780C3B', '#c5812a')
# pal_surv_dmet <- c('#8ac7ad', '#b74233')

dft_surv_dmet %<>%
  remove_trunc_gte_event(
    trunc_var = 'dmet_cpt_rep_yrs',
    event_var = 'tt_os_dmet_yrs'
  )

surv_obj_os_dmet_lt_adj <- with(
  dft_surv_dmet,
  Surv(
    time = dmet_cpt_rep_yrs,
    time2 = tt_os_dmet_yrs,
    event = os_dx_status # While it says dx, death is the same no matter the index time.
  )
)

# Get the median surv in months.  The rest is nonsense:
survfit(surv_obj_os_dmet_lt_adj ~ 1) %>%
  summary %>%
  `$`(., 'table') %>%
  multiply_by(., 12)

surv_obj_os_dmet_no_lt_adj <- with(
  dft_surv_dmet,
  Surv(
    time = tt_os_dmet_yrs,
    event = os_dx_status # While it says dx, death is the same no matter the index time.
  )
)

dft_km_no_lt_adj <- survfit(
  surv_obj_os_dmet_no_lt_adj ~ 1,
  data = dft_surv_dmet
) %>%
  broom::tidy(.) %>%
  add_row(time = 0, estimate = 1)


gg_os_dmet <- plot_one_survfit(
  dat = dft_surv_dmet,
  surv_form = surv_obj_os_dmet_lt_adj ~ 1,
  plot_title = "OS from metastasis",
  plot_subtitle = glue(
    "<span style = 'color:{pal_surv_dmet[1]};'>Adjusted</span> and
    <span style = 'color:{pal_surv_dmet[2]};'>Unadjusted</span> for left truncation (delayed entry)"
  ),
  x_exp = 0.1,
  force_color = pal_surv_dmet[1]
) +
  geom_step(
    data = dft_km_no_lt_adj,
    inherit.aes = F,
    aes(x = time, y = estimate),
    color = pal_surv_dmet[2]
  )


readr::write_rds(
  x = gg_os_dmet,
  file = here(surv_desc_fp, 'os_dmet_by_adjustment.rds')
)


# Test the independence of truncation and event times.
# There's not much we can do with this information, but it's a curiousity that
#   someone may be interested in.

# We'll just do all the tests and save them in case anyone asks.
# Rediculous practice statistically but we also don't plan to do anything with this.
dft_trunc_ind_test <- tribble(
  ~lab,
  ~dat,
  ~v_trunc,
  ~v_event,
  ~v_event_ind,
  "From diagnosis",
  dft_surv_dx,
  "dx_cpt_rep_yrs",
  "tt_os_dx_yrs",
  "os_dx_status",
  "From metastasis (all)",
  dft_surv_dmet,
  "dmet_cpt_rep_yrs",
  "tt_os_dmet_yrs",
  "os_dx_status",
  "From metastasis (metastatic at dx)",
  filter(dft_surv_dmet, stage_dx_iv %in% "Stage IV" & ca_dmets_yn %in% "Yes"),
  "dmet_cpt_rep_yrs",
  "tt_os_dmet_yrs",
  "os_dx_status",
) %>%
  slice(rep(1:n(), times = 3)) %>%
  mutate(method = rep(c("MB", "IPW1", "IPW2"), each = n() / 3))

dft_trunc_ind_test %<>%
  mutate(
    res = purrr::pmap(
      .l = list(
        dat = dat,
        var_trunc = v_trunc,
        var_event_time = v_event,
        var_event_ind = v_event_ind,
        method = method
      ),
      .f = test_kendall_tau
    )
  ) %>%
  select(lab, res) %>%
  unnest(res)

# IPW2 uses only the uncensored observations.  I don't really get why, so I'm confortable ignoring that it's different from the others.

readr::write_rds(
  x = dft_trunc_ind_test,
  file = here('data', 'survival', 'trunc_test.rds')
)

# Just a curiousity:  Unconditional.
# cor(
#   dft_surv_dmet$tt_os_dmet_yrs,
#   dft_surv_dmet$dmet_cpt_rep_yrs,
#   method = "kendall"
# )
