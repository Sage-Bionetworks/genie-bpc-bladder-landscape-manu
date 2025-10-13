# Description: Takes the oncotree codes from the CPT dataset and adds them
#   into the mutation, CNA and fusion files.  Also reshapes and filters the
#   CNA code so it runs in a reasonable amount of time.

library(purrr)
library(here)
library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

dft_mut <- readr::read_tsv(
  here('data-raw', 'genomic', 'data_mutations_extended.txt')
)
dft_cna <- readr::read_tsv(
  file = here("data-raw", "genomic", "data_CNA.txt")
)
dft_fus <- readr::read_tsv(
  file = here('data-raw', 'genomic', 'data_fusions.txt')
)
dft_cpt <- readr::read_rds(
  file = here('data', 'cohort', 'cpt.rds')
)

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
dft_otc <- dft_cpt %>%
  select(
    cpt_genie_sample_id,
    ONCOTREE_CODE = cpt_oncotree_code # required name for the annotator
  )
dft_mut <- left_join(
  dft_mut,
  dft_otc,
  by = c(Tumor_Sample_Barcode = "cpt_genie_sample_id"),
  relationship = "many-to-one"
)
dft_cna_long_selected <- left_join(
  dft_cna_long_selected,
  dft_otc,
  by = c(Tumor_Sample_Barcode = "cpt_genie_sample_id"),
  relationship = "many-to-one"
)
dft_fus <- left_join(
  dft_fus,
  dft_otc,
  by = c(Tumor_Sample_Barcode = "cpt_genie_sample_id"),
  relationship = "many-to-one"
)

# Spits out a message for each dataframe about the number of rows removed
#  and returns the data with the rows removed.
genomic_row_removed_helper <- function(dat) {
  dat_name <- deparse(substitute(dat))
  dat_nrow_pre <- nrow(dat)
  dat %<>% filter(!is.na(ONCOTREE_CODE))

  dat_nrow_post <- nrow(dat)
  nrow_diff <- dat_nrow_pre - dat_nrow_post
  nrow_diff_pct <- nrow_diff / dat_nrow_pre
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
  file = here('data', 'genomic', 'mut_ready_to_annotate.txt'),
  na = ''
)
readr::write_tsv(
  x = dft_cna_long_selected,
  file = here('data', 'genomic', 'cna_ready_to_annotate.txt'),
  na = ''
)
readr::write_tsv(
  x = dft_fus,
  file = here('data', 'genomic', 'fus_ready_to_annotate.txt'),
  na = ''
)
