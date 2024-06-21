# Create the co-occurence plots.  This is now a ton of code.
library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source) # also loads lots of packages.

out_dir <- here('data', 'genomic', 'gene_corr')


# Load in the needed data - probably a bit overkill here but it's copied.
read_wrap_clin <- function(p) {
  read_rds(file = here("data", 'cohort', p))
}

dft_pt <- read_wrap_clin("pt.rds")
dft_ca_ind <- read_wrap_clin("ca_ind.rds")
dft_cpt <- read_wrap_clin("cpt_aug.rds")

read_wrap_geno <- function(p) {
  read_rds(file = here("data", 'genomic', p))
}

dft_gp_all <- read_wrap_geno('gene_panel_all.rds')
dft_gp_sum <- read_wrap_geno('gene_panel_sum.rds')
dft_gp_by_gene <- read_wrap_geno('gene_panel_by_gene.rds')
dft_onco_impact <- read_wrap_geno('oncokb_impact.rds')
dft_alt <- read_wrap_geno('alterations.rds')
dft_gt_any_alt <- read_wrap_geno(
  'gene_test_any_alt.rds'
)



dft_inst_freq <- dft_gt_any_alt %>%
  select(
    sample_id,
    hugo,
    tested,
    any_alt,
    any_alt_onco
  ) %>%
  left_join(
    .,
    select(
      dft_cpt, 
      sample_id = cpt_genie_sample_id,
      institution
    ),
    by = "sample_id",
    relationship = 'many-to-one'
  )

dft_inst_freq %<>%
  group_by(institution, hugo) %>%
  summarize(
    n_tested = sum(tested, na.rm = T),
    n_alt = sum(any_alt, na.rm = T),
    n_alt_onco = sum(any_alt_onco, na.rm = T),
    
    prop_alt = n_alt/n_tested,
    prop_alt_onco = n_alt_onco/n_tested,
    .groups = "drop"
  ) 

dft_inst_freq_all <- dft_inst_freq %>%
  group_by(hugo) %>%
  summarize(
    across(
      .cols = c(n_alt, n_alt_onco, n_tested),
      .fns = sum
    ),
    .groups = 'drop'
  ) %>%
  mutate(
    prop_alt = n_alt/n_tested,
    prop_alt_onco = n_alt_onco/n_tested
  ) %>%
  mutate(institution = "All")





#########################################
# Build the whole-cohort co-occur plots #
#########################################

vec_co_occur_genes <- dft_inst_freq_all %>% 
  mutate(prop_tested = n_tested/max(n_tested)) %>%
  filter(prop_alt_onco > 0.1, prop_tested > 0.85) %>%
  pull(hugo)

dft_top_gene_bin <- dft_alt %>% 
  filter(hugo %in% vec_co_occur_genes) %>%
  filter(oncogenic %in% c("Likely Oncogenic", "Oncogenic")) %>%
  group_by(sample_id, hugo) %>%
  summarize(exists = n() >= 1, .groups = "drop") %>%
  pivot_wider(
    names_from = hugo,
    values_from = exists
  ) 

# Add in rows where nothing at all was found:
dft_top_gene_bin <- dft_cpt %>%
  select(sample_id = cpt_genie_sample_id) %>%
  left_join(dft_top_gene_bin, by = "sample_id")

dft_top_gene_bin %<>%
  mutate(
    across(
      .cols = -sample_id,
      .fns = (function(x) {
        x <- as.integer(x)
        x <- if_else(is.na(x), 0L, x)
      })
    )
  ) %>%
  # just ordering:
  select(sample_id, any_of(vec_co_occur_genes), everything())


dft_gene_assoc <- test_fisher_co_occur(
  dat = dft_top_gene_bin,
  ignore_cols = c("sample_id"),
  top = Inf, # already handled this above.
  alpha = 0.05
)

dft_gene_assoc %<>%
  mutate(
    assoc_lab = glue("{ct_11} / {ct_10+ct_01}"),
    p_value_adj = p.adjust(p.value, method = "BH"),
    cont_table_lab = cont_table_lab_help(
      ct_11 = ct_11,
      ct_01 = ct_01,
      ct_10 = ct_10,
      ct_00 = ct_00
    )
  )

gg_gene_assoc_sens <- plot_binary_association(
  dat = dft_gene_assoc,
  x_var = "var1_lab",
  y_var = "var2_lab",
  show_p_sig = T,
  label_var = NULL,
  pval_var = "p_value_adj"
)  + 
  theme(
    axis.text.x.top = element_text(angle = -20, hjust = 1)
  )

readr::write_rds(
  gg_gene_assoc_sens,
  # "All samples" here refers to tested or not.  Implicit assumption that 
  #    untested is negative.
  file = here('data', 'genomic', 'gene_corr', 'gg_mat_fisher_all_samples.rds')
)


# Sensitivity analysis to the above:  Only people who had all these genes tested in their panels.
# we want to force the same group of genes for the last test since
#   this is a sensitivity analysis.  This gets those in their axis
#. order:
vec_genes_in_co_occur_plot <- dft_gene_assoc$var1_lab %>%
  levels(.) %>%
  str_replace_all(
    ., "[:space:].*", ""
  )

vec_panels_with_all_top <- dft_gp_all %>%
  filter(
    hugo %in% vec_genes_in_co_occur_plot
  ) %>%
  select(cpt_seq_assay_id, hugo, tested) %>% 
  group_by(cpt_seq_assay_id) %>%
  summarize(
    has_all_top = n() >= length(vec_genes_in_co_occur_plot) # those with tested = F are implicitly missing currently.
  ) %>%
  filter(has_all_top) %>%
  pull(cpt_seq_assay_id) %>%
  as.character(.)

dft_alt_full_top_tested <- dft_alt %>%
  left_join(
    ., 
    select(dft_cpt, sample_id = cpt_genie_sample_id, cpt_seq_assay_id),
    by = 'sample_id'
  ) %>%
  filter(cpt_seq_assay_id %in% vec_panels_with_all_top)

# Now we copy-paste the same code from above:
dft_top_gene_bin_main <- dft_alt_full_top_tested %>% 
  filter(hugo %in% vec_genes_in_co_occur_plot) %>%
  filter(oncogenic %in% c("Likely Oncogenic", "Oncogenic")) %>%
  group_by(sample_id, hugo) %>%
  summarize(exists = n() >= 1, .groups = "drop") %>%
  pivot_wider(
    names_from = hugo,
    values_from = exists
  ) 

# Add in rows where nothing at all was found:
dft_top_gene_bin_main <- dft_cpt %>%
  filter(cpt_seq_assay_id %in% vec_panels_with_all_top) %>%
  select(sample_id = cpt_genie_sample_id) %>%
  left_join(dft_top_gene_bin, by = "sample_id")

dft_top_gene_bin_main %<>%
  mutate(
    across(
      .cols = -any_of('sample_id'),
      .fns = (function(x) {
        x <- as.integer(x)
        x <- if_else(is.na(x), 0L, x)
      })
    )
  ) %>%
  # just ordering:
  select(sample_id, vec_genes_in_co_occur_plot)

dft_gene_assoc_main <- test_fisher_co_occur(
  dat = dft_top_gene_bin_main,
  ignore_cols = c("sample_id"),
  top = Inf,
  alpha = 0.05,
  axis_order = vec_genes_in_co_occur_plot
)



dft_gene_assoc_main %<>%
  mutate(
    assoc_lab = glue("{ct_11} / {ct_10+ct_01}"),
    p_value_adj = p.adjust(p.value, method = "BH"),
    cont_table_lab = cont_table_lab_help(
      ct_11 = ct_11,
      ct_01 = ct_01,
      ct_10 = ct_10,
      ct_00 = ct_00
    )
  )

# This is now the "main" analysis.
gg_gene_assoc_main <- plot_binary_association(
  dat = dft_gene_assoc_main,
  x_var = "var1_lab",
  y_var = "var2_lab",
  show_p_sig = T,
  label_var = NULL,
  pval_var = "p_value_adj"
)  + 
  theme(
    axis.text.x.top = element_text(angle = 45, hjust = 0),
    plot.margin = unit(c(0.25, 1, 0.25, 0.25), "cm")
  )

readr::write_rds(
  gg_gene_assoc_main,
  file = here(out_dir, "gg_mat_fisher_tested_for_all_genes.rds")
)
  

# get the number of people in each cell from both plots for the sake of text explanation:
help_cell_table_total <- function(dat) {
  n_s <- dat %>%
    mutate(cell_n = ct_11 + ct_10 + ct_01 + ct_00) %>% 
    pull(cell_n)
  
  if (var(n_s, na.rm = T) > 0) {
    cli_abort("Problem:  Cells have different counts")
  } 
  return(unique(n_s))
}

n_cell_gene_assoc_sens <- help_cell_table_total(dft_gene_assoc)
n_cell_gene_assoc_main <- help_cell_table_total(dft_gene_assoc_main)



#################################
# Create plot comparing the two #
#################################

dft_gene_assoc_compare <-
  bind_rows(
    (dft_gene_assoc %>%
       filter(p_value_adj < 0.05) %>%
       select(var1, var2)),
    (dft_gene_assoc_main %>%
       filter(p_value_adj < 0.05) %>%
       select(var1, var2))
  ) %>%
  distinct(.)

dft_gene_assoc_compare <- dft_gene_assoc_main %>%
  select(
    var1, 
    var2, 
    main_or = odds_ratio,
    main_pval= p_value_adj
  ) %>%
  left_join(
    dft_gene_assoc_compare, .,
    by = c('var1', 'var2')
  )

dft_gene_assoc_compare <- dft_gene_assoc %>%
  select(
    var1, 
    var2, 
    sensitivity_or = odds_ratio,
    sensitivity_pval= p_value_adj
  ) %>%
  left_join(
    dft_gene_assoc_compare,
    .,
    by = c('var1', 'var2')
  )

dft_gene_assoc_compare %<>%
  pivot_longer(
    cols = matches("main|sensitivity")
  ) %>%
  separate(name, sep = "_", into = c("analysis", "quantity")) %>%
  pivot_wider(
    names_from = 'quantity',
    values_from = 'value'
  ) %>%
  mutate(pval_sig = if_else(pval < 0.05, T, F)) %>%
  arrange(desc(var1), desc(var2)) %>%
  mutate(lab = fct_inorder(paste0(var1, ":", var2)))

gg_gene_assoc_compare <- ggplot(
  dft_gene_assoc_compare,
  aes(y = lab, x = or, color = analysis, shape = pval_sig)
) + 
  geom_vline(xintercept = 1, linetype = "12") + 
  geom_point() + 
  scale_x_continuous(
    expand = expansion(mult = 0.01, add = 0),
    trans = 'log10',
    minor_breaks = log_tick_helper()
  ) +
  guides(x = guide_axis(minor.ticks = T)) +
  scale_color_vibrant() + 
  labs(
    title = "Estimate comparison for gene association analyses",
    subtitle = "Triangles = significant P value, Circles = not", 
    x = "Odds ratio (log10 scale)"
  ) + 
  guides(shape = "none") + 
  theme_bw() + 
  theme(
    axis.text.y = element_text(hjust = 0),
    axis.title.y = element_blank(),
    legend.position = "bottom",
    plot.title.position = "plot"
  )

readr::write_rds(
  gg_gene_assoc_compare,
  # "Single matrix" because we're about to split into metastatic and primary,
  #    which fans of math will know is two matrices.
  here('data', 'genomic', 'gene_corr', 'gg_sens_compare_single_matrix.rds')