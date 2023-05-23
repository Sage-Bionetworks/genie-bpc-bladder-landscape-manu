# Description: Grabs the raw data from Synapse and stores it in the data-raw folder.
# Author: Alex Paynter

# The Synapse folder containing the clinical data files.
synid_clin_data <- "syn28495599"  # "syn50612196"

library(cli)
library(synapser)
library(purrr)
library(dplyr)
library(here)
library(stringr)
library(magrittr)

synLogin()

# create directories for data and data-raw
dir.create(here("data"), showWarnings = F)
dir.create(here("data-raw"), showWarnings = F)

df_clin_children <- synGetChildren(synid_clin_data) %>%
  as.list %>%
  purrr::map_dfr(.x = .,
                 .f = as_tibble)

if (any(stringr::str_detect(df_clin_children$name, ".csv^"))) {
  warning("Non-CSV files unexpectedly contained in {synid_clin_data}.")
}

syn_store_in_dataraw <- function(sid) {
  synGet(entity = sid, downloadLocation = here("data-raw"))
}

purrr::walk(.x = df_clin_children$id, 
            .f = syn_store_in_dataraw)
