# ca_ind_dat is the index cancer diagnosis dataset.
make_dmet_status_block <- function(ca_ind_dat) {
  dmet_stat_levs <- c("No dmet noted", "Distant Metastasis")

  ci_sub <- ca_ind_dat %>%
    # just to be clear about variables we're using:
    select(
      record_id,
      ca_seq,
      stage_dx_iv,
      ca_dmets_yn,
      # dx_to_dmets_yrs, # not ever used, need helper.
      tt_os_dx_yrs
    )

  dmet_types <- get_dmet_time(ca_ind_dat, annotate_type = T)

  dmet_at_onset <- dmet_types %>%
    filter(.met_type %in% "stage_iv_with_immediate_dmet") %>%
    select(record_id, ca_seq) %>%
    left_join(., ci_sub, by = c("record_id", "ca_seq")) %>%
    mutate(
      dmet_status = dmet_stat_levs[2],
      dx_block_start = 0,
      dx_block_end = tt_os_dx_yrs
    )

  dmet_in_fu <- dmet_types %>%
    filter(!(.met_type %in% "stage_iv_with_immediate_dmet")) %>%
    # dx_dmet_yrs my derived version of the flawed dx_to_dmets_yrs.  See get_dmet_time() for construction.
    select(record_id, ca_seq, dx_dmet_yrs) %>%
    left_join(., ci_sub, by = c("record_id", "ca_seq")) %>%
    group_by(record_id, ca_seq) %>%
    # need two blocks for each of these people:
    slice(rep(1:n(), each = 2)) %>%
    mutate(
      dmet_status = if_else(
        row_number() %in% 1,
        dmet_stat_levs[1],
        dmet_stat_levs[2]
      ),
      dx_block_start = if_else(row_number() %in% 1, 0, dx_dmet_yrs),
      dx_block_end = if_else(row_number() %in% 1, dx_dmet_yrs, tt_os_dx_yrs)
    ) %>%
    ungroup(.)

  dmet_never <- anti_join(
    ci_sub,
    dmet_types,
    by = c("record_id", "ca_seq")
  ) %>%
    mutate(
      dmet_status = dmet_stat_levs[1],
      dx_block_start = 0,
      dx_block_end = tt_os_dx_yrs
    )

  dmet_all <- bind_rows(
    dmet_never,
    dmet_in_fu,
    dmet_at_onset
  ) %>%
    # keep tt_os_dx_yrs for sorting the plot.
    select(
      record_id,
      ca_seq,
      dmet_status,
      contains("dx_block"),
      tt_os_dx_yrs
    ) %>%
    mutate(dmet_status = factor(dmet_status, dmet_stat_levs))

  return(dmet_all)
}
