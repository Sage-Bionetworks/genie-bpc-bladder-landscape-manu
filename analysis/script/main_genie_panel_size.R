# this may be stored elsewhere but I couldn't find it.
# Just a simple count of the number of genes in each panel.

hugo_index <- readr::read_rds(
  here('data', 'genomic', 'main_genie', 'sample_hugo_index.rds')
)
samp_size <- hugo_index %>%
  filter(tested) %>%
  group_by(sample_id, patient_id, seq_assay_id, sample_type) %>%
  summarize(
    n_genes = length(unique(hugo)),
    .groups = 'drop'
  )
readr::write_rds(
  samp_size,
  here('data', 'genomic', 'main_genie', 'panel_size_by_sample.rds')
)
