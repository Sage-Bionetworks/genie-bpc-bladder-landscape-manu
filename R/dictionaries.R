# This file defines some dataframes - a rare exception to the rule that everything
#   in /R is a function.
dict_timing_vars <- tribble(
  ~raw,
  ~plot_long,
  "tt_first_reg_dx_yrs",
  "First regimen (from dx)",
  "tt_seq_dx_yrs",
  "First NGS (from dx)",
  "tt_seq_first_reg_yrs",
  "First NGS (from first regimen)",
  "tt_dmet_dx_yrs",
  "Distant Metastasis (from dx)",
  "tt_os_dx_yrs",
  "Followup, OS (from dx)",
) %>%
  mutate(default = plot_long)


dict_timing_vars_dmet <- tribble(
  ~raw,
  ~plot_long,
  "tt_seq_dmet_yrs",
  "First NGS (from dmet)",
  "tt_reg_post_dmet_yrs",
  "Post-dmet regimen (from dmet)",
  "tt_pfs_dmet_yrs",
  "Progression, PFS (from dmet)",
  "tt_os_dmet_yrs",
  "Followup, OS (from dmet)",
) %>%
  mutate(default = plot_long)
