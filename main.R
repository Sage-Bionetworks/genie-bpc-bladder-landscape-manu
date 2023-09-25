# Description: Top level workflow for the project.  Can be converted to a 
#   cleaner workflow later on.
# Author: Alex Paynter

library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)



############
# Clinical #
############
source(here('analysis', 'script', 'folder_setup.R'))
source(here('analysis', 'script', 'get_raw_data.R'))
# Create cohort data
# Create drug data
# Create dmet data
# Need to fill this in still.



###########
# Genomic #
###########
source(here('analysis', 'script', 'reshape_cna.R'))
source(here('analysis', 'script', 'prepare_data_for_oncokb_annotate.R'))
# Run annotate_oncokb.sh - see instructions there in comments.
source(here('analysis', 'script', 'create_gene_panel_dat.R'))
source(here('analysis', 'script', 'process_oncokb_output.R'))