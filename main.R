# Description: Top level workflow for the project.  Can be converted to a
#   cleaner workflow later on.
# Author: Alex Paynter

library(purrr)
library(here)
library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

# Manual downloads:
# - data-raw/manual/icd_topography.xlsx from https://www.ncri.ie/html/icdo3sites

############
# Clinical #
############
source(here('analysis', 'script', 'folder_setup.R'))
# get_raw_data.R takes several minutes, due to the need to get main GENIE data.
source(here('analysis', 'script', 'get_raw_data.R'))
source(here('analysis', 'script', 'create_cohort_data.R'))
source(here('analysis', 'script', 'create_characteristics_table.R'))
source(here('analysis', 'script', 'create_characteristics_table_main_genie.R'))
source(here('analysis', 'script', 'create_drug_dat.R'))
source(here('analysis', 'script', 'create_dmet_data.R'))
source(here('analysis', 'script', 'estimate_reg_year.R')) # resampling takes ~15s
source(here('analysis', 'script', 'create_neoadj_dat.R'))
source(here('analysis', 'script', 'derive_lines_of_therapy.R'))
source(here('analysis', 'script', 'derive_custom_time_to_met.R'))
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
cli_alert_warning(
  c(
    "The oncoprint needs to be manually created (sadly)\n",
    "follow the steps in \n/data/genomic/oncoprint/directions.md"
  )
)
# need to set up quarto render still.
rmarkdown::render(
  input = here('analysis', 'report', 'genie-bpc-bladder-genomic.qmd'),
  output_file = '02-genie-bpc-bladder-genomic.html',
  output_dir = here('output', 'report')
)


# Main genie genomic reprise
source(here(
  'analysis',
  'script',
  'prepare_data_for_oncokb_annotate_main_genie.R'
))
# Run annotate_oncokb_main_genie.sh - see instructions there in comments.
source(here('analysis', 'script', 'process_oncokb_output_main_genie.R'))
source(here('analysis', 'script', 'main_genie_panel_size.R'))

# DDR report:
source(here('analysis', 'script', 'oncoprint_ddr.R'))
source(here('analysis', 'script', 'ddr_definition_comparison.R'))
source(here('analysis', 'script', 'ddr_outcome_model_prep.R'))
source(here('analysis', 'script', 'ddr_outcome_model_run.R'))
source(here('analysis', 'script', 'ddr_outcome_panel_size_plots.R'))

# need to set up quarto render still.
quarto_render_copy <- function(
    qmd_path,
    qmd_file_stub,
    out_path,
    out_file = NULL
) {
  out_file <- out_file %||% paste0(qmd_file_stub, '.html')
  quarto::quarto_render(
    input = here(qmd_path, paste0(qmd_file_stub, '.qmd'))
  )
  fs::file_copy(
    path = here(qmd_path, paste0(qmd_file, '.html'))
    new_path = here(out_path, out_file),
    overwrite = TRUE
  )
}

quarto_render_copy(
  qmd_path = here('analysis', 'report'),
  qmd_file = here('genie-bpc_bladder_ddr'),
  out_path = here('output', 'report')
)
  
  
quarto::quarto_render(
  input = here('analysis', 'report', 'genie-bpc-bladder-ddr.qmd')
)
fs::file_copy(
  path = 'analysis/report/genie-bpc-bladder-ddr.html',
  new_path = 'output/report/02-genie-bpc-bladder-ddr.html',
  overwrite = TRUE
)
rmarkdown::render(
  input = here('analysis', 'report', 'genie-bpc-bladder-ddr.qmd'),
  output_file = '02-genie-bpc-bladder-ddr.html',
  output_dir = here('output', 'report')
)


############
# Survival #
############

source(here('analysis', 'script', 'basic_survival_descriptive.R'))
source(here('analysis', 'script', 'surv_ddr_data_all_1L.R'))
source(here('analysis', 'script', 'surv_ddr_onco_plat_univar.R'))
# Currently Figure 3 does not rely on anything but the univariate model.  Probably a concerning thing.
source(here('analysis', 'script', 'manu-fig3.R'))
source(here('analysis', 'script', 'manu-fig3-table-overlays.R'))
# The "interact" surv script is moved to deprecated and the model one needs to be updated to use the 1L data.

# Deleted the 2L script for now.  It can be created again by just changing the 1L script.

source(here('analysis', 'script', 'create_2025_asco_gu_abstract_figures.R'))
source(here('analysis', 'script', 'create_2025_asco_gu_abstract_figures_lot.R'))
source(here('analysis', 'script', 'create_2025_asco_gu_abstract_demo.R'))

source(here('analysis', 'script', 'create_fig_1.R'))


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
