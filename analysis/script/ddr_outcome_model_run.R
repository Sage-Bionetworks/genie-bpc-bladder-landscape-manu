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
      n_genes +
        de_novo_met +
        dob_ca_dx_yrs +
        dx_path_proc_cpt_yrs +
        upper_tract +
        birth_year +
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

ddr_out_all_mod <- bind_rows(
  ddr_out_mod_main_slr,
  ddr_out_mod_main_mlr
) %>%
  mutate(model = fct_inorder(model))

ddr_out_all_mod %<>% filter(!(term %in% "(Intercept)"))

write_rds(
  ddr_out_all_mod,
  here(dir_out, 'ddr_outcome_main_all_model_results_main.rds')
)

# just for plotting:
ddr_out_all_mod %<>% filter(!str_detect(term, "institution"))


gg_ddr_out_mod_compare <- forest_mod_natural_scale(
  dat = mutate(ddr_out_all_mod, model = fct_rev(model))
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
  gg_ddr_out_mod_compare,
  here(dir_out, 'gg_ddr_out_mod_compare_main.rds')
)
