# Description:  This script creates a figure showing the effect of adjusting
#   for truncation.  The point is a lunch and learn presentation at AACR.

library(purrr); library(fs); library(here)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

dft_ca_ind <- readr::read_rds(
  here('data', 'cohort', 'ca_ind.rds')
)

dft_cpt <- readr::read_rds(
  here('data', 'cohort', 'cpt.rds')
)

dft_first_cpt <- get_first_cpt(
  ca_ind_dat = dft_ca_ind,
  cpt_dat = dft_cpt
)

dft_surv_s4 <- dft_ca_ind %>%
  filter(stage_dx_iv %in% "Stage IV") %>%
  select(
    record_id, ca_seq,
    os_dx_status, tt_os_dx_yrs
  )

dft_surv_s4 %<>%
  left_join(., dft_first_cpt, by = c("record_id", "ca_seq"))

dft_surv_s4 %<>% filter(dx_cpt_rep_yrs < tt_os_dx_yrs)

surv_pal <- c('#bb5566', '#004488', '#ddaa33')
surv_pal_2 <- c("#4477aa", "#ee6677", "#228833", "#ccbb44", "#66ccee", "#aa3377", "#bbbbbb")

surv_obj_os_no_adjust <- with(
  dft_surv_s4,
  Surv(
    time = tt_os_dx_yrs,
    event = os_dx_status
  )
)

surv_obj_os <- with(
  dft_surv_s4,
  Surv(
    time = dx_cpt_rep_yrs,
    time2 = tt_os_dx_yrs,
    event = os_dx_status
  )
)

surv_obj_os_tm <- with(
  dft_surv_s4,
  trSurvfit(
    trun = dx_cpt_rep_yrs,
    obs = tt_os_dx_yrs,
    delta = os_dx_status,
    tFun = 'linear' # default
  )
)
dft_surv_os_tm <- surv_obj_os_tm$surv %>% 
  as_tibble %>%
  rename_all(tolower)

# We'll want the fit for this so we can plot it too:
dft_surv_os_no_adjust <- survfit(
  surv_obj_os_no_adjust ~ 1, data = dft_surv_s4
) %>%
  broom::tidy(.) %>%
  add_row(time = 0, estimate = 1)

dft_surv_os_ind_trunc <- survfit(
  surv_obj_os ~ 1, data = dft_surv_s4
) %>%
  broom::tidy(.) %>%
  add_row(time = 0, estimate = 1)





vec_surv_subtitle <- glue(
  "<span style = 'color:{surv_pal[1]};'>Not adjusted </span>and 
  <span style = 'color:{surv_pal[2]};'>Adjusted </span>for delayed entry (adjustment assumes independent truncation)"
)


gg_os <- survfit2(surv_obj_os ~ 1, data = dft_surv_s4) %>%
  ggsurvfit(color = surv_pal[2]) +
  geom_step(data = dft_surv_os_no_adjust,
            inherit.aes = F,
            aes(x = time, y = estimate),
            color = surv_pal[1]) + 
  # So tempting but it'll be confusing:
  # geom_step(data = dft_surv_os_tm,
  #           inherit.aes = F,
  #           aes(x = time, y = trsurv),
  #           color = surv_pal[3]) + 
  add_risktable(
    risktable_stats = c(
      "n.risk",
      "cum.censor",
      "cum.event"
    ),
    hjust = 0,
    # times = seq(2.5/2, 100, by = 2.5/2),
    size = 3.5 # default
  ) +
  scale_y_continuous(
    expand = c(0,0), 
    label = scales::label_percent()
  ) + 
  scale_x_continuous(
    name = "Years from cancer diagnosis",
    expand = c(0,0),
    breaks = seq(0, 100, by = 2.5)
  ) + 
  coord_cartesian(
    xlim = c(0, NA),
    ylim = c(0,1)
  ) + 
  labs(
    title = "Overall survival from diagnosis with Stage IV bladder cancer",
    subtitle = vec_surv_subtitle
  ) + 
  theme(
    axis.title.y = element_blank(),
    plot.title.position = "plot",
    plot.subtitle = element_markdown(
      margin = margin(c(0,0,15,0))
    ),
    # prevents the axis tick label clipping:
    plot.margin=unit(c(.2,.2,.2,.2),"cm")
  )

ggsave(
  filename = here('output', 'fig', 'lunch_n_learn_surv.jpeg'),
  plot = gg_os,
  dpi = 300,
  height = 5,
  width = 7
)
