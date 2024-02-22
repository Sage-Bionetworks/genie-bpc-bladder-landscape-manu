library(purrr); library(fs); library(here);
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

read_wrap <- function(p) {
  read_rds(
    file = here("data", 'cohort', paste0(p, '.rds'))
  )
}

dft_pt <- read_wrap("pt")
dft_med_onc <- read_wrap("med_onc")

dft_med_onc %>% pull(md_ecog) %>% unique %>% sort
dft_med_onc %>% pull(md_karnof) %>% unique %>% sort

dft_comp_counts <- dft_med_onc %>%
  mutate(
    ecog_complete = case_when(
      md_ecog %in% "Not documented in note" ~ 0,
      is.na(md_ecog) ~ 0,
      T ~ 1
    ),
    karnof_complete = case_when(
      md_karnof %in% "Not documented in the note" ~ 0,
      is.na(md_karnof) ~ 0,
      T ~ 1
    )
  ) %>%
  group_by(record_id) %>%
  summarize(
    n_complete_ecog = sum(ecog_complete),
    n_complete_karnof = sum(karnof_complete),
    .groups = 'drop'
  ) 

dft_comp_counts <- left_join(
  select(dft_pt, record_id),
  dft_comp_counts,
  by = 'record_id',
  relationship = 'one-to-one'
) 

dft_comp_counts %<>%
  mutate(
    across(
      .cols = contains("n_complete"),
      .fns = \(x) if_else(is.na(x), 0, x)
    )
  ) %>%
  mutate(
    n_either = n_complete_ecog + n_complete_karnof
  )

dft_comp_counts %>% summary


dft_comp_counts %>%
  summarize(
    n_karn_complete = sum(n_complete_karnof > 0),
    n_ecog_complete = sum(n_complete_ecog > 0),
    either_gt_zero = sum(n_complete_karnof > 0 | n_complete_ecog > 0)
  )
  
