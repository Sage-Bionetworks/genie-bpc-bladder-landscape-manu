make_grand_timing_df <- function(
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
      # will grab dmets separately
    ) %>%
    full_join(
      dft_timing, ., 
      by = c("record_id", "ca_seq")
    )
  
  dft_timing <- pt_dat %>% 
    select(record_id, birth_year) %>%
    full_join(
      dft_timing, .,
      by = c("record_id")
    )
  
  dft_timing <- get_dmet_time(ca_ind_dat) %>%
    # renaming variable to avoid downstream issues.
    select(record_id, ca_seq, tt_dmet_dx_yrs = dx_dmet_yrs) %>%
    full_join(
      dft_timing, .,
      by = c("record_id", "ca_seq")
    )

  dft_timing <- get_first_regimen(
    ca_ind_dat = ca_ind_dat,
    reg_dat = reg_dat
  ) %>%
    full_join(
      dft_timing, .,
      by = c("record_id", "ca_seq")
    )
  
  dft_timing %<>%
    # deviating from given names here for consistent formatting.
    rename(
      tt_dx_dob_yrs = dob_ca_dx_yrs,
      tt_seq_dx_yrs = dx_cpt_rep_yrs,
      tt_first_reg_dx_yrs  = dx_reg_start_int_yrs,
    ) %>%
    mutate(
      # Create age variables:
      age_dx = tt_dx_dob_yrs,
      age_first_reg = age_dx + tt_first_reg_dx_yrs,
      age_seq = age_dx + tt_seq_dx_yrs,
      age_dmet = age_dx + tt_dmet_dx_yrs,
      age_death = age_dx + tt_os_dx_yrs,
    ) %>%
    select(
      record_id, ca_seq,
      tt_dx_dob_yrs,
      tt_seq_dx_yrs,
      tt_first_reg_dx_yrs,
      tt_dmet_dx_yrs,
      tt_os_dx_yrs,
      
      age_dx,
      age_seq,
      age_first_reg,
      age_dmet,
      age_death
    )
  
  return(dft_timing)
}