# Description: Top level workflow for the project.  Can be converted to a 
#   cleaner workflow later on.
# Author: Alex Paynter

library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

# Manual downloads:
# - data-raw/manual/icd_topography.xlsx from https://www.ncri.ie/html/icdo3sites




############
# Clinical #
############
source(here('analysis', 'script', 'folder_setup.R'))
source(here('analysis', 'script', 'get_raw_data.R'))
source(here('analysis', 'script', 'create_cohort_data.R'))
source(here('analysis', 'script', 'create_characteristics_table.R'))
source(here('analysis', 'script', 'create_drug_dat.R'))
source(here('analysis', 'script', 'create_dmet_data.R'))
source(here('analysis', 'script', 'estimate_reg_year.R'))
source(here('analysis', 'script', 'create_neoadj_dat.R'))
source(here('analysis', 'script', 'derive_lines_of_therapy.R'))
rmarkdown::render(
  input = here('analysis', 'report', 'genie-bpc-bladder-clinical.rmd'),
  output_file = '01-genie-bpc-bladder-clinical.html',
  output_dir = here('output', 'report')
)
# Not used directly here, but could be later on:
source(here('analysis', 'script', 'create_time_invariant_model_factors.R'))
source(here('analysis', 'script', 'staging_clin_path.R'))




###########
# Genomic #
###########
source(here('analysis', 'script', 'prepare_data_for_oncokb_annotate.R'))
# Run annotate_oncokb.sh - see instructions there in comments.
source(here('analysis', 'script', 'create_gene_panel_dat.R'))
source(here('analysis', 'script', 'process_oncokb_output.R'))
source(here('analysis', 'script', 'add_tmb_to_cpt.R'))
source(here('analysis', 'script', 'plot_co_occur.R'))
# need to set up quarto render still.
rmarkdown::render(
  input = here('analysis', 'report', 'genie-bpc-bladder-genomic.qmd'),
  output_file = '02-genie-bpc-bladder-genomic.html',
  output_dir = here('output', 'report')
)
#  






############
# Survival #
############
source(here('analysis', 'script', 'basic_survival_descriptive.R'))
source(here('analysis', 'script', 'surv_platinum_gem_first_line.R'))
source(here('analysis', 'script', 'surv_platinum_gem_first_line_model.R'))
source(here('analysis', 'script', 'surv_ercc3_platinum.R'))
source(here('analysis', 'script', 'surv_ercc2_platinum.R'))
source(here('analysis', 'script', 'surv_first_line_immunotherapy.R'))
source(here('analysis', 'script', 'surv_ddr_onco_plat.R'))
source(here('analysis', 'script', 'surv_ddr_onco_plat_model.R'))
source(here('analysis', 'script', 'surv_ddr_neoadj.R'))
source(here('analysis', 'script', 'surv_second_line.R'))
# Render the qmd for survival.



#############################
# Met classification report #
#############################
# I cannot get this to render using quarto_render(), frustrating.
# Do:
# 1. render analysis/report/genie-bpc-met-class.qmd.
# 2. Run the following copy command:
fs::file_move(
  path = here('analysis', 'report', 'genie-bpc-bladder-met-class.html'),
  new_path = here('output', 'report', '99-genie-bpc-bladder-met-class.html')
)


