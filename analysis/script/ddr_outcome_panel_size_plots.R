library(purrr)
library(here)
library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)


dir_in <- here('data', 'genomic', 'ddr_def_compare', 'ddr_as_outcome')
dir_out <- here('data', 'genomic', 'ddr_def_compare', 'ddr_as_outcome')
ddr_outcome <- readr::read_rds(
  here(dir_in, 'ddr_outcome_mod_ready.rds')
)

ddr_outcome_main <- readr::read_rds(
  here(dir_in, 'ddr_outcome_mod_ready_main.rds')
)

inst_helper <- function(dat) {
  rtn <- dat %>%
    mutate(
      inst = str_extract(sample_id, "-[A-z]+-"),
      inst = str_sub(inst, 2, -2)
    )

  inst_lev <- c('DFCI', 'MSK', 'UHN', 'VICC', 'Other')

  rtn %>%
    mutate(
      inst_f = factor(inst, levels = inst_lev),
      inst_f = fct_explicit_na(inst_f, "Other")
    )
}


ddr_outcome_main %<>% inst_helper(.)
ddr_outcome %<>% inst_helper(.)

ddr_all <- bind_rows(
  (ddr_outcome_main %>%
    select(sample_id, ddr, seq_assay_id, panel_genes_100, inst_f) %>%
    mutate(study = "Main GENIE")),
  (ddr_outcome %>%
    select(sample_id, ddr, seq_assay_id, panel_genes_100, inst_f) %>%
    mutate(study = "BPC"))
) %>%
  mutate(study = fct_inorder(study))

ddr_all_aggr <- ddr_all %>%
  group_by(study, seq_assay_id, panel_genes_100, inst_f) %>%
  summarize(cases = n(), prop_ddr = mean(ddr), .groups = 'drop')


gg <- ggplot(
  ddr_all_aggr,
  aes(x = panel_genes_100, y = prop_ddr)
) +
  geom_point(aes(size = cases, color = inst_f, text = seq_assay_id)) +
  geom_smooth(
    aes(weight = cases),
    color = 'gray20',
    span = 1.5,
    se = F
  ) +
  geom_line(
    stat = 'smooth',
    method = 'lm',
    color = 'gray20',
    alpha = 0.5,
    linetype = '21',
    span = 1.5,
    se = F
  ) +
  theme_bw() +
  facet_wrap(vars(study), nrow = 1) +
  labs(
    x = "Genes tested in panel",
    y = "Percent DDR+"
  ) +
  scale_x_continuous(
    expand = expansion(add = c(0, 1), mult = 0),
    limits = c(0, NA),
    breaks = seq(0, 10, by = 2),
    labels = seq(0, 10 * 100, by = 2 * 100)
  ) +
  scale_y_continuous(
    expand = expansion(add = c(0.01, 0.05), mult = 0),
    limits = c(0, NA),
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_color_bmj()

# Plotly script is broken, which is OK since it's not permanent

readr::write_rds(
  gg,
  here(dir_out, 'pltly_ddr_panel_size.rds')
)
