library(purrr)
library(here)
library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

dir_out <- here('data', 'genomic', 'ddr_def_compare', 'ddr_as_outcome')

ddr_outcome <- readr::read_rds(
  here(
    'data',
    'genomic',
    'ddr_def_compare',
    'ddr_as_outcome',
    'ddr_outcome_mod_ready.rds'
  )
)

ddr_outcome %<>% select(-c(sample_id, seq_assay_id))
# Multiple imputation model first.

blocks <- construct.blocks(
  formulas = list(md_ecog_imp_num ~ .)
)
pm <- make.predictorMatrix(
  data = dplyr::select(
    ddr_outcome,
    -c(
      # taking away the outcome & id for predicting md_ecog_imp_num
      ddr
    )
  ),
  blocks = blocks
)

# If you treat ECOG scores as an ordinal variable, which is probably
#   more correct, you get 'polr' as the default method.  I this seems
#   fine for imputation, but it doesn't make a ton of sense to me for
#   the model (3 terms).

mids_surv <- mice(
  data = ddr_outcome,
  maxit = 5,
  m = 15,
  seed = 2341,
  predictorMatrix = pm,
  blocks = blocks,
  print = FALSE,
  # tended to be pretty similar to predictive mean matching:
  method = 'midastouch' # makes sense to me, gives similar results.
)

gg_imp <- plot_gg_strip(mids_surv, var = "md_ecog_imp_num", pt_size = 1) +
  labs(y = "Last observed ECOG")

readr::write_rds(
  x = gg_imp,
  file = here(dir_out, 'gg_imputations_ddr_outcome.rds')
)


# A bit painful here:  You can't use the data argument, or the "~." syntax.
#   If you do, then it will use the actual data without imputations.
#   I think there would be a way to autogen the formula and feed it in,
#   but it's only ~10 names so I'm going to type them out.
mira_surv <- with(
  mids_surv,
  glm(
    formula = ddr ~
      panel_genes_100 +
        de_novo_met +
        age_dx_5 +
        dx_path_proc_cpt_yrs +
        upper_tract +
        birth_year_5 +
        female +
        md_ecog_imp_num +
        met_sample +
        mibc_at_cpt +
        met_at_cpt +
        institution_DFCI +
        institution_UHN +
        institution_VICC +
        race_Asian +
        # removed the term for black - see note a few lines down.
        # race_Black +
        race_race_other_unk,
    family = 'binomial'
  )
)

ddr_out_mod_mi <- pool(mira_surv) %>%
  summary(conf.int = T) %>%
  select(term, estimate, conf.low = `2.5 %`, conf.high = `97.5 %`, p.value)
ddr_out_mod_mi %<>% mutate(model = "Multiple LR, imputed")

# Skipping the coefficient for Black people because it's totally colinear.
# Obviously a very important thing to note though - 0% of black people had a DDR alteration.
# ddr_outcome %>% filter(race_Black %in% 1)
ddr_out_mod_slr <- purrr::map_dfr(
  .x = colnames(select(ddr_outcome, -c(ddr, race_Black))),
  .f = \(x) {
    univariate_logistic_help(dat = ddr_outcome, x = x, y = 'ddr')
  }
) %>%
  mutate(model = "Simple LR")

ddr_out_all_mod <- bind_rows(
  ddr_out_mod_slr,
  ddr_out_mod_mi
) %>%
  mutate(model = fct_inorder(model))

ddr_out_all_mod %<>% filter(!(term %in% "(Intercept)"))

write_rds(
  ddr_out_all_mod,
  here(dir_out, 'ddr_outcome_all_model_results.rds')
)

# just for plotting:
ddr_out_all_mod %<>% filter(!str_detect(term, "institution"))


gg_ddr_out_mod_compare <- forest_mod_natural_scale(
  dat = mutate(ddr_out_all_mod, model = fct_rev(model))
) +
  labs(
    x = "Logistic regression coefficent<br>log(OR), original variable scale",
    y = NULL,
    title = "Associations with DDR alteration"
  ) +
  theme(
    plot.title.position = 'plot',
    axis.title = element_markdown()
  )

write_rds(
  gg_ddr_out_mod_compare,
  here(dir_out, 'gg_ddr_out_mod_compare.rds')
)


# Now do the same for main GENIE:

ddr_outcome_main <- readr::read_rds(
  here(
    'data',
    'genomic',
    'ddr_def_compare',
    'ddr_as_outcome',
    'ddr_outcome_mod_ready_main.rds'
  )
)

ddr_outcome_main %<>% select(-c(patient_id, sample_id, seq_assay_id))

ddr_out_mod_main_mlr <- glm(
  data = ddr_outcome_main,
  family = binomial,
  formula = ddr ~ .
)

ddr_out_mod_main_mlr <- ddr_out_mod_main_mlr %>%
  broom::tidy(., conf.int = T)

ddr_out_mod_main_mlr %<>% mutate(model = "Multiple LR")

ddr_out_mod_main_slr <- purrr::map_dfr(
  .x = colnames(select(ddr_outcome_main, -c(ddr))),
  .f = \(x) {
    univariate_logistic_help(dat = ddr_outcome_main, x = x, y = 'ddr')
  }
) %>%
  mutate(model = "Simple LR")

ddr_out_all_mod_main <- bind_rows(
  ddr_out_mod_main_slr,
  ddr_out_mod_main_mlr
) %>%
  mutate(model = fct_inorder(model))

ddr_out_all_mod_main %<>% filter(!(term %in% "(Intercept)"))

write_rds(
  ddr_out_all_mod_main,
  here(dir_out, 'ddr_outcome_main_all_model_results_main.rds')
)

# just for plotting:
ddr_out_all_mod_main %<>% filter(!str_detect(term, "institution"))


gg_ddr_out_mod_compare_main <- forest_mod_natural_scale(
  dat = mutate(ddr_out_all_mod_main, model = fct_rev(model))
) +
  labs(
    x = "Logistic regression coefficent<br>log(OR), original variable scale",
    y = NULL,
    title = "Associations with DDR alteration (main GENIE)"
  ) +
  theme(
    plot.title.position = 'plot',
    axis.title = element_markdown()
  )

# plotly::ggplotly(gg_ddr_out_mod_compare)

write_rds(
  gg_ddr_out_mod_compare_main,
  here(dir_out, 'gg_ddr_out_mod_compare_main.rds')
)


ddr_out_both <- bind_rows(
  mutate(ddr_out_all_mod_main, data = "Main GENIE"),
  mutate(ddr_out_all_mod, data = "BPC")
)

ddr_out_both %<>%
  # for simplicity we'll remove the imputed label (and put it in a legend)
  mutate(
    model = fct_collapse(
      model,
      `Multiple LR` = c("Multiple LR", "Multiple LR, imputed")
    )
  ) %>%
  arrange(
    model,
    data
  ) %>%
  mutate(
    multi_strat = fct_inorder(paste0(data, " - ", model))
  )

# So many weird things it's hard to automate this - we'll go manual
var_pretty_map <- tribble(
  ~raw,
  ~short,
  # uses the markdown format for plotting.
  "panel_genes_100",
  "**Panel size (/100 genes)**",

  "met_sampleTRUE",
  "Metastatic sample",

  'femaleTRUE',
  "Female",

  'upper_tractTRUE',
  "Upper tract UC",

  'race_race_other_unk',
  "Race = other/unknown",

  'race_Black',
  "Race = Black", # keeping even though BPC has no term.

  'race_Asian',
  "Race = Asian",

  "met_at_cptTRUE",
  "Metastatic @NGS",

  "de_novo_metTRUE",
  "de novo metastatic",

  "mibc_at_cptTRUE",
  "MIBC @NGS",

  "md_ecog_imp_num",
  "ECOG @NGS (imputed)",

  "dx_path_proc_cpt_yrs",
  "**Yrs from dx to path**",

  "birth_year_5",
  "**Birth year (/5 yrs)**",

  "age_seq_5",
  "**Age at seq (/5 yrs)**",

  "age_dx_5",
  "**Age at dx (/5 yrs)**",
)

ddr_out_both <- left_join(
  ddr_out_both,
  select(var_pretty_map, term = raw, term_lab = short),
  by = 'term'
) %>%
  mutate(
    term_lab = fct_rev(factor(term_lab, levels = var_pretty_map$short))
  )

ddr_out_both %<>%
  # for simplicity we'll remove the imputed label (and put it in a legend)
  mutate(
    model = fct_collapse(
      model,
      `Multiple LR` = c("Multiple LR", "Multiple LR, imputed")
    )
  ) %>%
  arrange(
    model,
    data
  ) %>%
  mutate(
    multi_strat = fct_inorder(paste0(data, " - ", model))
  )


gg_ddr_out_both_jumble <- forest_mod_natural_scale(
  ddr_out_both,
  model_var = "multi_strat",
  term_var = 'term_lab',
  pal = paste0('#', c('8fa6cc', 'd998c1', '2458ac', 'b20070'))
) +
  labs(
    x = "Coefficient (log(OR))",
    y = NULL,
    title = "Associations with DDR alteration",
    subtitle = "Bolded terms are continuous variables"
  ) +
  theme(
    plot.title.position = 'plot',
    axis.title = element_markdown()
  )

gg_ddr_out_both_panels <- forest_mod_natural_scale(
  ddr_out_both,
  model_var = "model",
  term_var = 'term_lab',
  pal = paste0('#', c('8fa6cc', 'b20070'))
) +
  labs(
    x = "Coefficient (log(OR))",
    y = NULL,
    title = "Associations with DDR alteration",
    subtitle = "Bolded terms are continuous variables"
  ) +
  theme(
    plot.title.position = 'plot',
    axis.title = element_markdown()
  ) +
  facet_wrap(vars(data))


gg_ddr_out_both_multiple_lr <- forest_mod_natural_scale(
  filter(ddr_out_both, model %in% "Multiple LR"),
  model_var = "data",
  term_var = 'term_lab',
  pal = paste0('#', c('2458ac', 'b20070'))
) +
  labs(
    x = "Coefficient [log(OR)]",
    y = NULL,
    title = "Associations with DDR alteration - multiple regression only",
    subtitle = "Bolded = continuous variables"
  ) +
  theme(
    plot.title.position = 'plot',
    axis.title = element_markdown()
  )

bundle_both <- list(
  ddr_out_both_data = list(ddr_out_both),
  gg_ddr_out_both_jumble = list(gg_ddr_out_both_jumble),
  gg_ddr_out_both_panels = list(gg_ddr_out_both_panels),
  gg_ddr_out_both_multiple_lr = list(gg_ddr_out_both_multiple_lr)
)

ggsave(
  plot = gg_ddr_out_both_panels,
  height = 4,
  width = 8,
  filename = here(
    'output',
    'fig',
    'aacr_summer_2025',
    'ddr_outcome_main_bpc.png'
  )
)

write_rds(
  bundle_both,
  here(dir_out, 'bundle_bpc_main_combined.rds')
)
