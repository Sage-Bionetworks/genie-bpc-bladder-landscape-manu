make_grand_timing_df_dmet <- function(
  pt_dat,
  ca_ind_dat,
  cpt_dat,
  reg_dat
) {
  dft_timing <- get_first_cpt(ca_ind_dat, cpt_dat)

  dft_timing <- ca_ind_dat %>%
    select(
      record_id,
      ca_seq,
      dob_ca_dx_yrs,
      tt_os_dx_yrs,
      tt_pfs_i_and_m_adv_yrs
    ) %>%
    full_join(
      dft_timing,
      .,
      by = c("record_id", "ca_seq")
    )

  dft_timing <- pt_dat %>%
    select(record_id, birth_year) %>%
    full_join(
      dft_timing,
      .,
      by = c("record_id")
    )

  dft_timing <- get_dmet_time(ca_ind_dat) %>%
    # renaming variable to avoid downstream issues.
    select(record_id, ca_seq, tt_dmet_dx_yrs = dx_dmet_yrs) %>%
    full_join(
      dft_timing,
      .,
      by = c("record_id", "ca_seq")
    )

  dft_timing <- get_first_regimen(
    ca_ind_dat = ca_ind_dat,
    reg_dat = reg_dat
  ) %>%
    full_join(
      dft_timing,
      .,
      by = c("record_id", "ca_seq")
    )

  dft_timing %<>%
    # deviating from given names here for consistent formatting.
    rename(
      tt_dx_dob_yrs = dob_ca_dx_yrs,
      tt_seq_dx_yrs = dx_cpt_rep_yrs,
      tt_first_reg_dx_yrs = dx_reg_start_int_yrs,
    ) %>%
    mutate(
      tt_dmet_dob_yrs = tt_dx_dob_yrs + tt_dmet_dx_yrs,
      tt_seq_dmet_yrs = tt_seq_dx_yrs - tt_dmet_dx_yrs,
      tt_reg_post_dmet_yrs = tt_first_reg_dx_yrs - tt_dmet_dx_yrs,
      # PFS is already relative to dmet, does not need to be subtracted off.
      tt_pfs_dmet_yrs = tt_pfs_i_and_m_adv_yrs,
      tt_os_dmet_yrs = tt_os_dx_yrs - tt_dmet_dx_yrs
    ) %>%
    select(
      record_id,
      ca_seq,
      tt_dmet_dob_yrs,
      tt_seq_dmet_yrs,
      tt_reg_post_dmet_yrs,
      tt_pfs_dmet_yrs,
      tt_os_dmet_yrs
    )

  chk_reg <- dft_timing %>%
    filter(tt_reg_post_dmet_yrs <= 0) %>%
    nrow

  if (chk_reg > 0) {
    cli::cli_abort(
      message = "There are regimens before the dmet event - function not designed for this use.  Please pass in a filtered version of the regimen data."
    )
  }

  return(dft_timing)
}
