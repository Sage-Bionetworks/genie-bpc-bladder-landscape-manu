





dft_inst_freq_all %>% 
  mutate(prop_tested = n_tested/max(n_tested)) %>%
  filter(prop_alt > 0.1, prop_tested > 0.85)

dft_top_gene_bin <- dft_alt %>% 
  # filter(hugo %in% vec_alt_most_common[1:20]) %>%
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
  select(sample_id, all_of(vec_alt_most_common), everything())


dft_gene_assoc <- test_fisher_co_occur(
  dat = dft_top_gene_bin,
  ignore_cols = c("sample_id"),
  top = 20,
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

gg_gene_assoc <- plot_binary_association(
  dat = dft_gene_assoc,
  x_var = "var1_lab",
  y_var = "var2_lab",
  show_p_sig = T,
  label_var = "cont_table_lab",
  pval_var = "p_value_adj"
)  + 
  theme(
    axis.text.x.top = element_text(angle = -20, hjust = 1)
  )













dft_inst_freq_all %>% 
  mutate(prop_tested = n_tested/max(n_tested)) %>%
  filter(prop_alt > 0.1, prop_tested > 0.85)

dft_top_gene_bin <- dft_alt %>% 
  # filter(hugo %in% vec_alt_most_common[1:20]) %>%
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
  select(sample_id, all_of(vec_alt_most_common), everything())


dft_gene_assoc <- test_fisher_co_occur(
  dat = dft_top_gene_bin,
  ignore_cols = c("sample_id"),
  top = 20,
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

gg_gene_assoc <- plot_binary_association(
  dat = dft_gene_assoc,
  x_var = "var1_lab",
  y_var = "var2_lab",
  show_p_sig = T,
  label_var = "cont_table_lab",
  pval_var = "p_value_adj"
)  + 
  theme(
    axis.text.x.top = element_text(angle = -20, hjust = 1)
  )












# Putting the code for the sensitivity here.  I want to make sure I get both.


# Sensitivity analysis to the above:  Only people who had all these genes tested in their panels.
# we want to force the same group of genes for the last test since
#   this is a sensitivity analysis.  This gets those in their axis
#. order:
vec_genes_in_co_occur_plot <- dft_gene_assoc$var1_lab %>%
  levels(.) %>%
  str_replace_all(
    ., "[:space:].*", ""
  )

vec_panels_with_all_top10 <- dft_gp_all %>%
  filter(
    hugo %in% vec_genes_in_co_occur_plot
  ) %>%
  select(cpt_seq_assay_id, hugo, tested) %>% 
  group_by(cpt_seq_assay_id) %>%
  summarize(
    has_all_top_10 = n() >= 10 # those with tested = F are implicitly missing currently.
  ) %>%
  filter(has_all_top_10) %>%
  pull(cpt_seq_assay_id) %>%
  as.character(.)

dft_alt_full_top_10_tested <- dft_alt %>%
  left_join(
    ., 
    select(dft_cpt, sample_id = cpt_genie_sample_id, cpt_seq_assay_id),
    by = 'sample_id'
  ) %>%
  filter(cpt_seq_assay_id %in% vec_panels_with_all_top10)

# Now we copy-paste the same code from above:
dft_top_gene_bin_sens <- dft_alt_full_top_10_tested %>% 
  filter(hugo %in% vec_genes_in_co_occur_plot) %>%
  group_by(sample_id, hugo) %>%
  summarize(exists = n() >= 1, .groups = "drop") %>%
  pivot_wider(
    names_from = hugo,
    values_from = exists
  ) 

# Add in rows where nothing at all was found:
dft_top_gene_bin_sens <- dft_cpt %>%
  filter(cpt_seq_assay_id %in% vec_panels_with_all_top10) %>%
  select(sample_id = cpt_genie_sample_id) %>%
  left_join(dft_top_gene_bin, by = "sample_id")

dft_top_gene_bin_sens %<>%
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
  select(sample_id, vec_genes_in_co_occur_plot)

dft_gene_assoc_sens <- test_fisher_co_occur(
  dat = dft_top_gene_bin_sens,
  ignore_cols = c("sample_id"),
  top = 10,
  alpha = 0.05,
  axis_order = vec_genes_in_co_occur_plot
)



dft_gene_assoc_sens %<>%
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
  dat = dft_gene_assoc_sens,
  x_var = "var1_lab",
  y_var = "var2_lab",
  show_p_sig = T,
  label_var = "cont_table_lab",
  pval_var = "p_value_adj"
)  + 
  theme(
    axis.text.x.top = element_text(angle = -20, hjust = 1)
  )



```



```{r}
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

n_cell_gene_assoc <- help_cell_table_total(dft_gene_assoc)
n_cell_gene_assoc_sens <- help_cell_table_total(dft_gene_assoc_sens)

n_cell_gene_assoc
n_cell_gene_assoc_sens