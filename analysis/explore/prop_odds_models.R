library(purrr)
library(here)
library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

dft_ddr_neo <- readr::read_rds(
  here('data', 'survival', 'ddr_neoadj', 'ddr_neoadj.rds')
)
dft_mod <- dft_ddr_neo %>%
  select(record_id, onco_ddr, clin_group_clust, path_group_clust) %>%
  pivot_longer(
    cols = contains("clust"),
    names_to = "t_type",
    values_to = "group_stage"
  ) %>%
  mutate(
    group_stage = forcats::fct_drop(group_stage),
    group_stage_num = as.numeric(group_stage),
    t = if_else(t_type %in% "clin_group_clust", 1, 2)
  )


library(MASS)
party <- factor(rep(c("Rep", "Dem"), c(407, 428)), levels = c("Rep", "Dem"))
rpi <- c(30, 46, 148, 84, 99) # cell counts
dpi <- c(80, 81, 171, 41, 55) # cell counts
ideology <- c(
  "Very Liberal",
  "Slightly Liberal",
  "Moderate",
  "Slightly Conservative",
  "Very Conservative"
)
pol.ideology <- factor(
  c(rep(ideology, rpi), rep(ideology, dpi)),
  levels = ideology
)
dat <- data.frame(party, pol.ideology)

# fit proportional odds model
library(MASS)
pom <- polr(pol.ideology ~ party, data = dat)


# Apply this simply to our case:

#install.packages('ordinal')
library(ordinal)
fm1 <- clm(rating ~ temp * contact, data = wine)
fm1 ## print method
summary(fm1)
broom::tidy(fm1, conf.int = T)

broom::tidy(fm1, conf.int = F) |> # checking my construction
  mutate(
    conf.low = estimate - std.error * qnorm(0.975),
    conf.high = estimate + std.error * qnorm(0.975)
  )


fm1 <- clm(rating ~ temp, data = wine, nominal = ~contact)

fm2 <- update(fm1, ~ . - temp:contact)
anova(fm1, fm2)
drop1(fm1, test = "Chi")
add1(fm1, ~ . + judge, test = "Chi")
fm2 <- step(fm1)
summary(fm2)


# This is actually worth using:
clmm(
  data = dft_mod,
  group_stage ~ t * onco_ddr + (1 | record_id)
) |>
  broom::tidy(conf.int = T) %>%
  # put the estimates and CIs on the odds ratio scale:
  mutate(
    across(
      .cols = c(estimate, conf.low, conf.high),
      .fns = \(z) exp(z)
    )
  )

# Just curious if I'm right about coarsening this being a bad idea:

dft_mod %<>%
  mutate(
    group_stage_coarse = forcats::fct_collapse(
      group_stage,
      III = c("IIIA(N0)", "IIIA(N1)", "IIIB")
    )
  )

clmm(
  data = dft_mod,
  group_stage_coarse ~ t * onco_ddr + (1 | record_id)
) |>
  broom::tidy(conf.int = T) %>%
  # put the estimates and CIs on the odds ratio scale:
  mutate(
    across(
      .cols = c(estimate, conf.low, conf.high),
      .fns = \(z) exp(z)
    )
  )
# Pretty similar.

# A second model for checking the sanity of this (plots too):

dft_mod %<>%
  mutate(
    group_stage_coarse = forcats::fct_collapse(
      group_stage,
      III = c("IIIA(N0)", "IIIA(N1)", "IIIB")
    )
  )

clmm(
  data = dft_mod,
  group_stage_coarse ~ t * onco_ddr + (1 | record_id)
) |>
  broom::tidy(conf.int = T) %>%
  # put the estimates and CIs on the odds ratio scale:
  mutate(
    across(
      .cols = c(estimate, conf.low, conf.high),
      .fns = \(z) exp(z)
    )
  )

dft_mod_2 <- dft_ddr_neo %>%
  dplyr::select(record_id, onco_ddr, group_stage_change_f) %>%
  mutate(
    group_stage_change_f = factor(
      group_stage_change_f,
      ordered = T,
      levels = levels(group_stage_change_f)
    )
  )

polr(data = dft_mod_2, group_stage_change_f ~ onco_ddr) %>%
  broom::tidy(., conf.int = T)
