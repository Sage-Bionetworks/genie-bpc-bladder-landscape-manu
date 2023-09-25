# Description: Grabs the raw data from Synapse and stores it in the data-raw folder.
# Author: Alex Paynter

library(purrr); library(fs); library(here)
purrr::walk(.x = fs::dir_ls('R'), .f = source)

# The Synapse folder containing the clinical data files.
synid_clin_data <- "syn28495599"
synid_cbio_data <- "syn26958249"

# genomic files to grab (panels are all grabbed based on file name):
geno_files_included <- c(
  "data_mutations_extended.txt",
  "data_CNA.txt",
  "data_fusions.txt",
  "data_cna_hg19.seg"
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

syn_store_in_dataraw <- function(sid) {
  synGet(entity = sid, downloadLocation = here("data-raw"))
}

purrr::walk(.x = df_clin_children$id, 
            .f = syn_store_in_dataraw)







# Get the genomic data from the "cBioPortal_files" directory.
df_geno_children <- synGetChildren(synid_cbio_data) %>%
  as.list %>%
  purrr::map_dfr(.x = .,
                 .f = as_tibble)

df_geno_children %<>%
  mutate(
    is_panel = str_detect(name, "^data_gene_panel_.*\\.txt$"),
    is_included = name %in% geno_files_included
  ) %>%
  filter(is_panel | is_included) 

syn_store_in_dataraw_geno <- function(sid) {
  synGet(entity = sid, downloadLocation = here("data-raw", "genomic"))
}

purrr::walk(.x = df_geno_children$id, 
            .f = syn_store_in_dataraw_geno)
