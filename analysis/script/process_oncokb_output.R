# Description:  Assess the impact of oncoKB filtering, save files with
#  only the variants which pass an oncoKB pass.

library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

dft_mut_onco <- readr::read_tsv(
  here('data', 'genomic', 'mut_onco.txt'),
  show_col_types = F
)
dft_cna_onco <- readr::read_tsv(
  here('data', 'genomic', 'cna_onco.txt'),
  show_col_types = F
)
dft_fus_onco <- readr::read_tsv(
  here('data', 'genomic', 'fus_onco.txt'),
  show_col_types = F
)

# For CNAs we focus only on the highly amplified cases, taking the decisions from the breast manuscript as a good starting point if nothing else.
dft_cna_raw <- readr::read_tsv(
  here('data', 'genomic', 'cna_ready_to_annotate.txt'),
  show_col_types = F
) 
dft_cna_onco <- bind_cols(
  # The sample ID that came back is garbage.
  select(dft_cna_onco, -SAMPLE_ID),
  select(dft_cna_raw,  SAMPLE_ID = Tumor_Sample_Barcode)
) 

anno_msg_help <- function(dat) {
  dat_name <- deparse(substitute(dat))
  
  tab <- tabyl(dat, ANNOTATED) %>%
    filter(ANNOTATED) # T/F so this grabs the annotated pct line.
  
  cli::cli_alert_info(
    text = glue("{dat_name}: Of the {tab$n} rows, {round(tab$percent*100,1)}% were annotated.")
  )
}

# Just some print outs for the analyst - hoping for near 100% in all:
anno_msg_help(dft_mut_onco)
anno_msg_help(dft_cna_onco)
anno_msg_help(dft_fus_onco)

onco_count_help <- function(dat, label) {
  tabyl(dat, ONCOGENIC) %>%
    mutate(type = label) %>%
    select(type, oncogenic = ONCOGENIC, n)
}

dft_onco_impact <- bind_rows(
  onco_count_help(dft_mut_onco, "Mutation"),
  onco_count_help(dft_cna_onco, "CNA"),
  onco_count_help(dft_fus_onco, "Fusion")
)

lev_onco <- c("Oncogenic", "Likely Oncogenic",
              "Likely Neutral",
              "Inconclusive", "Unknown")

dft_onco_impact %<>%
  mutate(
    oncogenic = factor(oncogenic, levels = lev_onco),
    type = fct_inorder(type)
  )

readr::write_rds(
  x = dft_onco_impact,
  file = here('data', 'genomic', 'oncokb_impact.rds')
)




# Addon: we want to calculate the variant allele frequency also.
# t_ref_count is the "read depth" (number of times a specific seqeunce was read)
#   supporting the reference sequence in this sample.
# t_alt_count is the read depth for the variant sequence in this sample.
# these both come from the BAM file, so even though we don't have the reads we 
# can still estimate VAF.
# A very small number of samples have this missing, which seems odd (but fine 
#   if it's really just two samples.

dft_mut_onco %<>%
  mutate(
    mut_vaf = t_alt_count / (t_alt_count + t_ref_count)
  )




# Create an alterations dataset - one row per alteration.

dft_mut_onco_alt <- dft_mut_onco %>%
  rename_all(tolower) %>%
  mutate(alt_type = "Mutation") %>%
  select(
    sample_id = tumor_sample_barcode,
    hugo = hugo_symbol,
    alt_type,
    hgvsc,
    hgvsp,
    chromosome,
    oncogenic,
    mutation_effect,
    highest_level,
    consequence,
    variant_classification,
    variant_type,
    mutation_status,
    mut_vaf
  )

dft_cna_onco_alt <- dft_cna_onco %>% 
  rename_all(tolower) %>%
  mutate(alt_type = "CNA") %>%
  select(
    sample_id,
    hugo = hugo_symbol,
    alt_type,
    cna_desc = alteration,
    oncogenic,
    highest_level,
    mutation_effect
  )

dft_fus_onco_alt <- dft_fus_onco %>% 
  rename_all(tolower) %>%
  mutate(alt_type = "Fusion") %>%
  select(
    sample_id = tumor_sample_barcode,
    hugo = hugo_symbol,
    alt_type,
    fusion_desc = fusion,
    oncogenic,
    mutation_effect,
    highest_level
  )

dft_alt <- bind_rows(
  dft_mut_onco_alt,
  dft_cna_onco_alt,
  dft_fus_onco_alt
)

# We need a unique key.  Currently a sample can have several alterations in the same #  hugo code and alteration type.  
# Initially I wanted to use the descriptions of the
#  alterations, such as the HGVSc code or the fusion description.  Because these
#  are sometimes missing I'm going to assign a number instead.  Sometimes even with 
#  a missing HGVSc&HGVSp code the variant can be annotated by oncokb - I don't 
#  want to remove those just to fit what would have been a beautiful data structure.
dft_alt %<>% 
  group_by(sample_id) %>%
  mutate(alt_seq = 1:n()) %>%
  ungroup(.) 

lev_alt_type <- c(
  "Mutation",
  "CNA",
  "Fusion"
)


# bit of cleanup
dft_alt %<>%
  mutate(
    alt_type = factor(alt_type, levels = lev_alt_type),
    oncogenic = factor(oncogenic, levels = lev_onco),
    mutation_effect = format_mutation_effect(mutation_effect),
    mut_eff_simple = format_mutation_effect_simple(mutation_effect),
    highest_level = format_highest_level(highest_level)
  ) %>%
  select(
    sample_id,
    alt_seq,
    hugo,
    alt_type,
    # features common to all alterations:
    oncogenic,
    mutation_effect,
    mut_eff_simple,
    highest_level,
    # descriptors for each alteration type:
    hgvsc,
    hgvsp,
    cna_desc,
    fusion_desc,
    # Any other data we pulled, like mutation stuff or whatever:
    everything()
  )

readr::write_rds(
  x = dft_alt,
  file = here('data', 'genomic', 'alterations.rds')
)


dft_gp_all <- readr::read_rds(
  here('data', 'genomic', 'gene_panel_all.rds')
)

gene_gte1_alt <- dft_alt %>% 
  pull(hugo) %>% unique %>% sort

gene_mut_not_in_panel <- setdiff(
  (dft_alt %>% filter(alt_type %in% "Mutation") %>%
     pull(hugo) %>% unique %>% sort),
  (dft_gp_all$hugo %>% unique %>% sort)
)

if (length(gene_mut_not_in_panel) > 0) {
  cli::cli_alert_danger(
    "Mutations exist which are in no panels: {paste(gene_mut_not_in_panel, collapse = ', ')}"
  )
}

readr::write_rds(
  gene_gte1_alt,
  file = here('data', 'genomic', 'gene_gte1_alt.rds')
)





# Assess the oncoKB impact on individual mutations

dft_cpt <- readr::read_rds(
  here('data', 'cohort', 'cpt.rds')
)
dft_gene_test <- dft_cpt %>%
  select(
    cpt_genie_sample_id, record_id, ca_seq, cpt_seq_assay_id, 
    contains("sample_type")
  )
dft_gp_all <- readr::read_rds(
  here('data', 'genomic', 'gene_panel_all.rds')
) %>%
  mutate(hugo = as.character(hugo))

dft_gene_test %<>%
  left_join(., dft_gp_all, by = "cpt_seq_assay_id",
            relationship = "many-to-many")


# Side venture
dft_alt_any <- dft_alt %>% 
  group_by(sample_id, hugo) %>%
  summarize(
    any_alt = T,
    any_alt_onco = (
      sum(
        oncogenic %in% c("Likely Oncogenic", "Oncogenic"),
        na.rm = T
      ) 
      >= 1),
    .groups = "drop"
  ) 

dft_gene_test_any_alt <- left_join(
  rename(dft_gene_test, sample_id = cpt_genie_sample_id),
  dft_alt_any,
  by = c("sample_id", "hugo")
) 

dft_gene_test_any_alt %<>%
  select(sample_id, record_id, ca_seq, cpt_seq_assay_id, sample_type,
         hugo, tested, any_alt, any_alt_onco) %>%
  mutate(
    any_alt = if_else(is.na(any_alt), F, any_alt),
    any_alt_onco = if_else(is.na(any_alt_onco), F, any_alt_onco)
  )

# put the dataset keys first:
dft_gene_test_any_alt %<>%
  select(sample_id, hugo, everything())

readr::write_rds(
  dft_gene_test_any_alt,
  file = here('data', 'genomic', 'gene_test_any_alt.rds')
)





dft_gene_test %<>%
  select(
    sample_id = cpt_genie_sample_id, 
    cpt_seq_assay_id, 
    hugo, 
    contains("tested")
  ) %>%
  pivot_longer(
    cols = contains("tested"),
    names_to = "test_type",
    values_to = "test_logical"
  ) %>%
  mutate(
    alt_type = case_when(
      test_type %in% "tested" ~ "Mutation",
      test_type %in% "tested_cna" ~ "CNA",
      test_type %in% "tested_fusion" ~ "Fusion"
    ),
    alt_type = factor(alt_type, levels = lev_alt_type)
  )

dft_gene_test %<>%
  filter(test_logical) %>%
  select(-c(test_type, test_logical))




dft_gene_test <- dft_alt %>% 
  # everything in this dataframe represents an alteration, so:
  mutate(altered = 1) %>%
  select(sample_id, hugo, alt_type, altered, oncogenic, mut_eff_simple) %>%
  left_join(
    dft_gene_test,
    .,
    by = c("sample_id", "hugo", "alt_type")
  )

# Not 100% sure this is needed, but it can't hurt:
readr::write_rds(
  x = dft_gene_test,
  file = here('data', 'genomic', 'gene_test_index.rds')
)


n_sample <- nrow(dft_cpt)
dft_prop_samples_gene_tested <- dft_gene_test %>%
  count(sample_id, hugo) %>% # n here is meaningless.
  count(hugo) %>%
  mutate(prop = n/n_sample) %>%
  arrange(desc(prop))

readr::write_rds(
  x = dft_prop_samples_gene_tested,
  file = here('data', 'genomic', 'gene_prop_samp_test.rds')
)



# Now we're ready to summarize:

dft_onco_imp <- dft_gene_test %>%
  group_by(hugo, alt_type) %>%
  summarize(
    n_tested = n(),
    n_altered = sum(altered, na.rm = T),
    n_oncogenic = sum(oncogenic %in% c("Likely Oncogenic", "Oncogenic"), na.rm = T),
    n_gain = sum(mut_eff_simple %in% "Gain", na.rm = T),
    n_loss = sum(mut_eff_simple %in% "Loss", na.rm = T),
    n_switch = sum(mut_eff_simple %in% "Switch", na.rm = T),
    .groups = "drop"
  )

readr::write_rds(
  x = dft_onco_imp,
  file = here('data', 'genomic', 'gene_counts.rds')
)

















