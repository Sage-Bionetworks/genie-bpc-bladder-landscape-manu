# Description: Takes the oncotree codes from the CPT dataset and adds them
#   into the mutation, CNA and fusion files.  Also reshapes and filters the
#   CNA code so it runs in a reasonable amount of time.

library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

dir_input <- here('data-raw', 'genomic', 'main_genie')
dir_output <- here('data', 'genomic', 'main_genie')

dft_mut <- fread(
  here(dir_input, 'data_mutations_extended.txt')
)
dft_cna <- fread(
  file = here(dir_input, "data_CNA.txt")
)
dft_fus <- fread(
  file = here(dir_input, 'data_fusions.txt')
)
clin_samp <- fread(
  file = here(dir_input, 'data_clinical_sample.txt'),
  skip = 4
)

dft_mut %<>% as_tibble(.)
dft_cna %<>% as_tibble(.)
dft_fus %<>% as_tibble(.)
clin_samp %<>% as_tibble(.) %>% rename_all(tolower)

# Pulled from the data guide for bladder:
onco_bladder <- c('BLAD', 'BLCA', 'BLSC', 'UTUC', 'SCBC')

samp_bladder <- clin_samp %>%
  filter(oncotree_code %in% onco_bladder) %>%
  pull(sample_id)

readr::write_rds(
  samp_bladder,
  here(dir_output, 'bladder_samples_mg.rds')
)

dft_mut %<>% filter(Tumor_Sample_Barcode %in% samp_bladder)
dft_cna %<>% select(Hugo_Symbol, any_of(samp_bladder))
dft_fus %<>% filter(Tumor_Sample_Barcode %in% samp_bladder)


# The CNA file needs to be reshaped to be fed into the annotator the way I know how:
dft_cna_long_selected <- dft_cna %>% 
  pivot_longer(
    cols = -Hugo_Symbol,
    names_to = "Tumor_Sample_Barcode", # just to match.
    values_to = "value"
  ) %>%
  # Trusting in the work of my collaborators 100% here:
  filter(!is.na(value) & value >= 2) 

# For each of our three file types we will add in the oncotree code:
dft_otc <- clin_samp %>%
  filter(sample_id %in% samp_bladder) %>%
  select(
    sample_id,
    ONCOTREE_CODE = oncotree_code
  )

dft_mut <- left_join(
  dft_mut, dft_otc, 
  by = c(Tumor_Sample_Barcode = "sample_id"),
  relationship = "many-to-one"
)
dft_cna_long_selected <- left_join(
  dft_cna_long_selected, dft_otc,
  by = c(Tumor_Sample_Barcode = "sample_id"),
  relationship = "many-to-one"
)
dft_fus <- left_join(
  dft_fus, dft_otc,
  by = c(Tumor_Sample_Barcode = "sample_id"),
  relationship = "many-to-one"
)

# Spits out a message for each dataframe about the number of rows removed
#  and returns the data with the rows removed.
genomic_row_removed_helper <- function(dat) {
  dat_name <- deparse(substitute(dat))
  dat_nrow_pre <- nrow(dat)
  dat %<>% filter(!is.na(ONCOTREE_CODE))
  
  dat_nrow_post <- nrow(dat)
  nrow_diff <- dat_nrow_pre-dat_nrow_post
  nrow_diff_pct <- nrow_diff/dat_nrow_pre
  cli::cli_alert_success(
    "Removed {nrow_diff} rows ({round(nrow_diff_pct,0)}%) from {dat_name} filtering down to only 'cohort' samples with a valid oncotree code."
  )
  
  return(dat)
}

dft_mut <- genomic_row_removed_helper(dft_mut)
dft_cna <- genomic_row_removed_helper(dft_cna_long_selected)
dft_fus <- genomic_row_removed_helper(dft_fus)


readr::write_tsv(
  x = dft_mut,
  file = here(dir_output, 'mut_ready_to_annotate.txt'),
  na = ''
)
readr::write_tsv(
  x = dft_cna_long_selected,
  file = here(dir_output, 'cna_ready_to_annotate.txt'),
  na = ''
)
readr::write_tsv(
  x = dft_fus,
  file = here(dir_output, 'fus_ready_to_annotate.txt'),
  na = ''
)
