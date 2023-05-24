output_synid <- "syn51536061" #2023-04-17-PrCa-landscape-paper-outputs

library(synapser)
library(magrittr)
library(here)

synLogin()
synapser::File(here("analysis", "reports", "bpc_bladder.html"),
               parent = output_synid) %>%
  synStore()
