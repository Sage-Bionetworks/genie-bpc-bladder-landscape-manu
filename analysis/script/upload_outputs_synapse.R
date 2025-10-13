output_synid <- "syn51536061"

library(synapser)
library(magrittr)
library(here)

synLogin()
synapser::File(
  here("analysis", "reports", "bpc_bladder.html"),
  parent = output_synid
) %>%
  synStore()
