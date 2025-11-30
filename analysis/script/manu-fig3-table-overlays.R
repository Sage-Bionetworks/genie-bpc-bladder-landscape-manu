library(purrr)
library(here)
library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source) # also load

manu_out_dir_fig3 <- here('output', 'manu', 'fig3')

met_plat_1l <- readr::read_rds(
  here('data', 'survival', 'ddr_onco_1L', 'met_ddr_surv_plat_1L.rds')
)

met_plat_1l %<>%
  remove_trunc_gte_event(
    trunc_var = 'reg_fcpt_yrs',
    event_var = 'tt_os_g_yrs'
  )

met_plat_1l %<>%
  mutate(reg_fcpt_yrs = ifelse(reg_fcpt_yrs < 0, 0, reg_fcpt_yrs))

met_plat_1l %<>%
  mutate(
    carbo = str_detect(regimen_drugs, "arboplatin"),
    cis = str_detect(regimen_drugs, 'isplatin')
  ) %>%
  filter(!(carbo & cis)) %>% # one case, not real (0 length)
  mutate(
    ddr_plat_comb = case_when(
      ddr_before_entry & carbo ~ "Carbo, DDR+",
      ddr_before_entry & cis ~ "Cis, DDR+",
      carbo ~ "Carbo, DDR-",
      cis ~ "Cis, DDR-"
    )
  ) %>%
  mutate(ddr_plat_comb = factor(ddr_plat_comb))


met_plat_1l_grouped <- bind_rows(
  (met_plat_1l %>% mutate(analysis_group = "any_plat")),
  (met_plat_1l %>%
    filter(str_detect(regimen_drugs, "Cisplatin")) %>%
    mutate(analysis_group = "cisplatin_based")),
  (met_plat_1l %>%
    filter(str_detect(regimen_drugs, "Carboplatin")) %>%
    mutate(analysis_group = "carboplatin_based"))
)

met_plat_1l_grouped %<>%
  nest(.by = analysis_group)


met_plat_1l_grouped %<>%
  mutate(
    cox_ddr = purrr::map(
      .x = data,
      .f = \(z) {
        coxph(
          data = z,
          with(
            z,
            Surv(
              time = reg_fcpt_yrs,
              time2 = tt_os_g_yrs,
              event = os_g_status
            )
          ) ~
            ddr_before_entry
        ) %>%
          broom::tidy(., conf.int = T, exponentiate = F) # will exp later.
      }
    )
  )


met_plat_1l_grouped %<>%
  mutate(
    median_surv = purrr::map(
      .x = data,
      .f = \(z) {
        survfit(
          data = z,
          with(
            z,
            Surv(
              time = reg_fcpt_yrs,
              time2 = tt_os_g_yrs,
              event = os_g_status
            )
          ) ~
            ddr_before_entry
        ) %>%
          summary %>%
          `$`(., 'table') %>%
          as_tibble(., rownames = 'group') %>%
          rename(lower = `0.95LCL`, upper = `0.95UCL`)
      }
    )
  )

# Removed a printout here.

fig3_cox_tab <- met_plat_1l_grouped %>%
  select(analysis_group, cox_ddr) %>%
  unnest(cox_ddr) %>%
  # put the estimates on the HR scale (not log HR):
  mutate(
    across(
      .cols = c(estimate, conf.low, conf.high),
      .fns = exp
    )
  ) %>%
  select(1, 3, 6:8) %>%
  mutate(
    ci_str = est_int_str(
      estimate,
      conf.low,
      conf.high,
      plus_prefix = F,
      est_digits = 2,
      na = "NA"
    )
  )


fig3_km_tab <- met_plat_1l_grouped %>%
  select(analysis_group, median_surv) %>%
  unnest(median_surv) %>%
  # put the estimates in months:
  mutate(
    across(
      .cols = c(median, lower, upper),
      .fns = \(z) z * 12.0148 # previously 12, makes almost no difference.
    )
  ) %>%
  # for now we just care about the medians:
  select(1, 2, n = records, median, lower, upper) %>%
  mutate(
    ci_str = est_int_str(
      median,
      lower,
      upper,
      plus_prefix = F,
      est_digits = 1,
      na = "NA"
    )
  )

# This function basically pulls the row data for one model and one group.
tab3a_help <- function(
  a_grp,
  grp,
  include_p = F
) {
  last_block <- (fig3_cox_tab %>%
    filter(analysis_group %in% a_grp) %>%
    pull(ci_str))
  if (include_p) {
    last_block <- paste0(
      last_block,
      '; ',
      p = format.pval(
        fig3_cox_tab %>%
          filter(analysis_group %in% a_grp) %>%
          pull(p.value),
        digits = 1,
        eps = 0.001
      )
    )
  }

  c(
    (fig3_km_tab %>%
      filter(
        analysis_group %in% a_grp,
        group %in% grp
      ) %>%
      pull(n)),
    (fig3_km_tab %>%
      filter(
        analysis_group %in% a_grp,
        group %in% grp
      ) %>%
      pull(ci_str)),
    last_block
  )
}

fig3a_table_inlay <- tibble(
  left_head = c(
    'n',
    'Median OS (95%)',
    'HR'
  ),
  DDRa = tab3a_help('ddr_before_entry=TRUE', a_grp = 'any_plat', include_p = T),
  `No DDRa` = tab3a_help(
    'ddr_before_entry=FALSE',
    a_grp = 'any_plat',
    include_p = T
  )
)

fig3a_table_inlay %>%
  rename(` ` = left_head) %>%
  flextable(.) %>%
  flextable::merge_at(i = 3, j = 2:3) %>%
  autofit(.) %>%
  save_as_docx(path = here(manu_out_dir_fig3, 'fig_3a_inlay_table.docx'))


fig3b_top <- tibble(
  left_head = c(
    'n',
    'Median OS (95%)',
    'HR'
  ),
  DDRa = tab3a_help('ddr_before_entry=TRUE', a_grp = 'cisplatin_based'),
  `No DDRa` = tab3a_help('ddr_before_entry=FALSE', a_grp = 'cisplatin_based')
) %>%
  mutate(left_grp_head = 'cisplatin') %>%
  relocate(left_grp_head, .before = left_head)

fig3b_bot <- tibble(
  left_head = c(
    'n',
    'Median OS (95%)',
    'HR'
  ),
  DDRa = tab3a_help('ddr_before_entry=TRUE', a_grp = 'carboplatin_based'),
  `No DDRa` = tab3a_help('ddr_before_entry=FALSE', a_grp = 'carboplatin_based')
) %>%
  mutate(left_grp_head = 'carboplatin') %>%
  relocate(left_grp_head, .before = left_head)

fig3b_table_inlay <- bind_rows(
  fig3b_top,
  fig3b_bot
)


fig3b_table_inlay %>%
  rename(`  ` = left_grp_head, ` ` = left_head) %>%
  flextable(.) %>%
  flextable::merge_v(j = 1) %>%
  flextable::valign(valign = "top") %>%
  flextable::merge_at(i = 3, j = 3:4) %>%
  flextable::merge_at(i = 6, j = 3:4) %>%
  autofit(.) %>%
  save_as_docx(path = here(manu_out_dir_fig3, 'fig_3b_inlay_table.docx'))
