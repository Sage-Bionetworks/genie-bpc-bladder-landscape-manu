# This file defines some dataframes - a rare exception to the rule that everything
#   in /R is a function.
dict_timing_vars <- tribble(
  ~raw, ~plot_long,
  "tt_seq_dx_yrs", "First NGS (from dx)",
  "tt_first_reg_dx_yrs", "First regimen (from dx)",
  "tt_dmet_dx_yrs", "Distant Metastasis (from dx)",
  "tt_os_dx_yrs", "Followup, OS (from dx)"
) %>%
  mutate(default = plot_long)
