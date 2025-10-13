# Description: Grab all the gene panel files and create data frames
#   to organize that information.

library(purrr)
library(here)
library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

vec_gene_panels <- fs::dir_ls(here('data-raw', 'genomic')) %>%
  str_filter(., 'data_gene_panel_.*')

gp_all <- purrr::map_dfr(
  .x = vec_gene_panels,
  .f = tidy_gene_panel
)

# for consistency with clinical data (and I like the name better)
gp_all %<>% rename(cpt_seq_assay_id = stable_id)


# Merge in the sample data.
dft_cpt <- readr::read_rds(
  here('data', 'cohort', 'cpt.rds')
)

gp_sum <- dft_cpt %>%
  select(
    cpt_genie_sample_id,
    record_id,
    cpt_number,
    sample_type,
    cpt_seq_assay_id
  ) %>%
  group_by(cpt_seq_assay_id) %>%
  summarize(
    n_pts = length(unique(record_id)),
    n_samples = n(),
    .groups = "drop"
  )

# Make sure we have a gene panel data file for each one found in the CPT data:
chk_gp_files <- setdiff(
  unique(gp_sum$cpt_seq_assay_id),
  unique(gp_all$cpt_seq_assay_id)
)
if (length(chk_gp_files) > 0) {
  cli_abort(
    glue(
      "Gene panels exist in CPT data which are not in the meta files:  {paste(chk_gp_files, collapse = ',')}"
    )
  )
}


# get the number of genes in each panel, merge that in to the summary dataframe.
gp_sum <- gp_all %>%
  group_by(cpt_seq_assay_id) %>%
  summarize(
    # n() would also be fine here:
    n_genes = length(unique(hugo)),
    .groups = 'drop'
  ) %>%
  left_join(
    gp_sum,
    .,
    by = "cpt_seq_assay_id"
  )
gp_sum %<>%
  arrange(desc(n_pts)) %>%
  mutate(
    cpt_seq_assay_id = forcats::fct_inorder(cpt_seq_assay_id)
  )

saveRDS(
  object = gp_sum,
  file = here('data', 'genomic', 'gene_panel_sum.rds')
)


# Metadata from the data guide (main GENIE data guide, table 4) which
#  I'm manually typing in:
# Specifically this is gene-level CNA and 'structural variants' are fusions.
dft_gp_data_guide <- tribble(
  ~panel,
  ~tested_cna,
  ~tested_fusion,
  "MSK-IMPACT341",
  1,
  1,
  "MSK-IMPACT410",
  1,
  1,
  "MSK-IMPACT468",
  1,
  1,
  "DFCI-ONCOPANEL-1",
  1,
  1,
  "DFCI-ONCOPANEL-2",
  1,
  1,
  "DFCI-ONCOPANEL-3",
  1,
  1,
  # Making a guess here since it doesn't show up explicitly:
  "DFCI-ONCOPANEL-3.1",
  1,
  1,
  "VICC-01-T7",
  1,
  1,
  "VICC-01-T5A",
  1,
  1,
  "VICC-01-SOLIDTUMOR",
  1,
  1, # bit of a guess here too.
  "UHN-48-V1",
  0,
  0,
  "UHN-50-V2",
  0,
  0,
  "UHN-555-V1",
  0,
  0,
  "UHN-555-BLADDER-V1",
  0,
  0,
  "UHN-OCA-V3",
  0,
  1
) %>%
  mutate(
    across(
      .cols = c(tested_cna, tested_fusion),
      .fns = as.logical
    )
  )

# Check that all rows in the gene data have a match here.
# Important with copy-pasted code between cohorts.
chk_data_guide_complete <- setdiff(
  unique(gp_all$cpt_seq_assay_id),
  unique(dft_gp_data_guide$panel)
)
if (length(chk_data_guide_complete)) {
  cli::cli_abort(
    glue(
      "Need to update dft_gp_data_guide to include all represented panels.  Not covered: {paste(chk_data_guide_complete, collapse = ',')}"
    )
  )
}

gp_all %<>%
  left_join(
    .,
    dft_gp_data_guide,
    by = c(cpt_seq_assay_id = "panel")
  )

gp_by_gene <- gp_all %>%
  group_by(hugo) %>%
  summarize(
    n_panels = length(unique(cpt_seq_assay_id)),
    .groups = "drop"
  )

# Find the proportion of samples which include testing for each gene.
gp_by_gene_samp_counts <- dft_cpt %>%
  select(cpt_genie_sample_id, cpt_seq_assay_id) %>%
  mutate(n_samples = n()) %>%
  left_join(
    .,
    select(gp_all, cpt_seq_assay_id, hugo),
    by = "cpt_seq_assay_id",
    relationship = "many-to-many"
  ) %>%
  group_by(hugo) %>%
  summarize(
    num_samp_tested = n(),
    prop_samp_tested = n() / first(n_samples),
    .groups = "drop"
  )

gp_by_gene <- left_join(gp_by_gene, gp_by_gene_samp_counts, by = "hugo")

gp_by_gene %<>%
  arrange(desc(n_panels), hugo) %>%
  mutate(
    hugo = forcats::fct_inorder(hugo)
  )

saveRDS(
  object = gp_by_gene,
  file = here('data', 'genomic', 'gene_panel_by_gene.rds')
)


# apply the factor levels above to the "all" data:
gp_all %<>%
  mutate(
    cpt_seq_assay_id = factor(
      cpt_seq_assay_id,
      levels = levels(gp_sum$cpt_seq_assay_id)
    ),
    hugo = factor(
      hugo,
      levels = levels(gp_by_gene$hugo)
    )
  )

saveRDS(
  object = gp_all,
  file = here('data', 'genomic', 'gene_panel_all.rds')
)
