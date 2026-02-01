# Description: Grabs the raw data from Synapse and stores it in the data-raw folder.
# Author: Alex Paynter

library(purrr)
library(fs)
library(here)
purrr::walk(.x = fs::dir_ls('R'), .f = source)

# The Synapse folder containing the clinical data files.
synid_clin_data <- "syn28495599"
synid_cbio_data <- "syn26958249"
# These have to be pulled in from main GENIE:
synid_assay_info <- 'syn22159815'
synid_bed_file <- 'syn9734427'
synid_bed_file_version <- 108


# genomic files to grab (panels are all grabbed based on file name):
geno_files_included <- c(
  "data_mutations_extended.txt",
  "data_CNA.txt",
  "data_fusions.txt",
  'data_sv.txt'
  # "data_cna_hg19.seg"
  # don't see a genomic information file
)

synLogin()


df_clin_children <- synGetChildren(synid_clin_data) %>%
  as.list %>%
  purrr::map_dfr(
    .x = .,
    .f = as_tibble
  )

if (any(stringr::str_detect(df_clin_children$name, ".csv^"))) {
  warning("Non-CSV files unexpectedly contained in {synid_clin_data}.")
}

syn_store_in_dataraw <- function(sid, loc = here('data-raw'), v = NULL) {
  # not sure how to do this with synGet, so we'll do a conditional for the version.
  if (is.null(v)) {
    synGet(
      entity = sid,
      downloadLocation = loc,
      ifcollision = 'overwrite.local'
    )
  } else {
    synGet(
      entity = sid,
      downloadLocation = loc,
      ifcollision = 'overwrite.local',
      version = v
    )
  }
}

purrr::walk(.x = df_clin_children$id, .f = syn_store_in_dataraw)


# Get the genomic data from the "cBioPortal_files" directory.
df_geno_children <- synGetChildren(synid_cbio_data) %>%
  as.list %>%
  purrr::map_dfr(.x = ., .f = as_tibble)

df_geno_children %<>%
  mutate(
    is_panel = str_detect(name, "^data_gene_panel_.*\\.txt$"),
    is_included = name %in% geno_files_included
  ) %>%
  filter(is_panel | is_included)

purrr::walk(.x = df_geno_children$id, .f = \(z) {
  syn_store_in_dataraw(z, loc = here("data-raw", "genomic"))
})


syn_store_in_dataraw(
  synid_assay_info,
  loc = here('data-raw', 'genomic')
)
# Insanely slow, probably due partially to the zillions of extra rows.
syn_store_in_dataraw(
  synid_bed_file,
  loc = here('data-raw', 'genomic'),
  v = synid_bed_file_version
)


# This now stays firmly on 18.2 consortium using hard codes.
# Using 18.2 because that's what we used for the pan-GENIE DDR analysis.
# Tired of version shifting and errors.
mg_clin_pt <- 'syn9734568.234' # notation means: synid.version
mg_clin_samp <- 'syn9734573.244'
mg_geno_info <- 'syn7444851.192'
mg_mut <- 'syn5571527.332'
mg_cna <- 'syn7213997.234'
mg_sv <- 'syn44810233.33'
mg_ai <- 'syn21614837.92'

purrr::walk(
  .x = c(mg_clin_pt, mg_clin_samp, mg_geno_info, mg_mut, mg_cna, mg_sv, mg_ai),
  .f = \(z) {
    syn_store_in_dataraw(
      z,
      loc = here("data-raw", "genomic", 'main_genie')
    )
  }
)
