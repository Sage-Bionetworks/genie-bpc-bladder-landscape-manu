# Description:  Assess the impact of oncoKB filtering, save files with
#  only the variants which pass an oncoKB pass.

library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

dir_input <- here('data', 'genomic', 'main_genie')
dir_output <- dir_input

dft_mut_onco <- readr::read_tsv(
  here(dir_input, 'mut_onco.txt'),
  show_col_types = F
)
dft_cna_onco <- readr::read_tsv(
  here(dir_input, 'cna_onco.txt'),
  show_col_types = F
)
dft_fus_onco <- readr::read_tsv(
  here(dir_input, 'fus_onco.txt'),
  show_col_types = F
)

# For CNAs we focus only on the highly amplified cases, taking the decisions from the breast manuscript as a good starting point if nothing else.
dft_cna_raw <- readr::read_tsv(
  here(dir_input, 'cna_ready_to_annotate.txt'),
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
  file = here(dir_output, 'oncokb_impact.rds')
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



# Also add in the DDR pathways
dft_pearl_pathways <- readr::read_rds(
  here('data', 'genomic', 'pearl_pathways.rds')
)

lev_path_pearl <- c(
  (dft_pearl_pathways$pathway %>% unique),
  "None"
)
dft_alt %<>% 
  left_join(
    .,
    select(dft_pearl_pathways, hugo = gene, pathway_ddr = pathway),
    by = "hugo",
    relationship = 'many-to-one'
  ) %>%
  mutate(
    pathway_ddr = if_else(is.na(pathway_ddr), "None", pathway_ddr),
    pathway_ddr = factor(pathway_ddr, levels = lev_path_pearl)
  )


readr::write_rds(
  x = dft_alt,
  file = here(dir_output, 'alterations.rds')
)

# Hoping that I don't need any of the rest of the script - see the NON main 
#   genie version for that if you need to pull from it.

# The one piece we will need is the dataframe listing who was tested for what genes.
# We will use the bed file to do that here.

samp_bladder <- read_rds(here(dir_input, 'bladder_samples_mg.rds'))
clin_samp_bladder <- readr::read_tsv(
  here('data-raw', 'genomic', 'main_genie', 'data_clinical_sample.txt'),
  skip = 4
) %>% rename_all(tolower) %>%
  filter(sample_id %in% samp_bladder)

panels_used <- clin_samp_bladder %>% pull(seq_assay_id) %>% unique(.)

mg_bed <- readr::read_tsv(
  here('data-raw', 'genomic', 'main_genie', 'genomic_information.txt')
) 

testing_df <- mg_bed %<>% 
  rename_all(tolower) %>%
  filter(includeinpanel %in% T) %>%
  filter(seq_assay_id %in% panels_used)

testing_df <- mg_bed %>%
  filter(!is.na(hugo_symbol)) %>%
  select(seq_assay_id, hugo_symbol) %>%
  distinct(.) %>%
  mutate(tested = T) %>%
  tidyr::complete(seq_assay_id, hugo_symbol, fill = list(tested = F))

testing_df %>% filter(seq_assay_id %in% "WAKE-CLINICAL-DX1", hugo_symbol %in% "TERT")

sample_hugo_index <- clin_samp_bladder %>% 
  select(sample_id, patient_id, seq_assay_id, sample_type)

sample_hugo_index <- 
  left_join(
    sample_hugo_index,
    testing_df,
    by = 'seq_assay_id',
    relationship = 'many-to-many'
  ) 

sample_hugo_index %<>% rename(hugo = hugo_symbol)

sample_hugo_index <- dft_alt %>% 
  group_by(sample_id, hugo) %>%
  summarize(
    any_alt = T,
    # for QC:
    any_mut = any(alt_type %in% "Mutation"),
    any_alt_onco = any(oncogenic %in% c('Oncogenic', 'Likely Oncogenic')),
    .groups = 'drop'
  ) %>%
  left_join(
    sample_hugo_index,
    .,
    by = c('sample_id', 'hugo')
  )

# QC curiousity:
# sample_hugo_index %>%
#   filter(!tested & any_mut) %>%
#   count(seq_assay_id, hugo, sort = T)
# The TERT examples are interesting.  The id says TERT and the hugo symbol
#   does not, meaning the site thought it was on TERT and the Sage remap
#   say it isn't.  This brings up an interesting question: Why do we only
#   remap the bed file regions and not the results?  We're releasing TERT
#   results even when the bed file indicates no TERT coverage - insanity.

# For now we will do a hideous hack and assume that any altered sample
#   was tested.  Which is true, but does not actually adjust correctly for
#   many samples which were not altered.
sample_hugo_index %<>%
  filter(any_mut) %>%
  select(seq_assay_id, hugo) %>%
  distinct(.) %>%
  mutate(panel_has_pos_test = T) %>%
  left_join(
    sample_hugo_index,
    .,
    by = c('seq_assay_id', 'hugo')
  )

n_row_sh <- nrow(sample_hugo_index)
n_hacked <- sample_hugo_index %>%
  filter(!tested & panel_has_pos_test) %>% nrow
sample_hugo_index %<>%
  mutate(tested = case_when(
    panel_has_pos_test ~ T,
    T ~ F
  ))
cli::cli_alert_warning(
  "{n_hacked} of the {n_row_sh} rows were changed to indicate testing based on a positive result in the panel for that hugo symbol (despite the bed file not indicating testing).  Not ideal but it prevents paradoxes.")

# This is a fairly large file to save, but it's the core of what we're doing.
readr::write_rds(
  sample_hugo_index,
  here(dir_output, 'sample_hugo_index.rds')
)
  

hugo_sum <- sample_hugo_index %>%
  group_by(hugo) %>%
  summarize(
    n_tested = n(),
    n_tested_chk = length(unique(sample_id)), # will remove.
    n_alt = length(unique(sample_id[any_alt])),
    n_alt_onco = length(unique(sample_id[any_alt_onco])),
  ) 

if (any(hugo_sum$n_tested != hugo_sum$n_tested_chk)) { 
  cli_abort("sample_hugo_index does not have unique sample/hugo rows.")
} else {
  hugo_sum %<>% select(-n_tested_chk)
}

hugo_sum %<>%
  mutate(
    prop_alt = n_alt / n_tested,
    prop_alt_onco = n_alt / n_tested
  )
    
readr::write_rds(
  hugo_sum,
  here(dir_output, 'hugo_summary.rds')
)





  













