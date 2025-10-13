#' @param dat_med_onc A medical oncology dataset from PRISSMM.
#' @param dat_cutoff A dataset of cutoff times (maxima) for relevant med onc notes.
#' @param var_cutoff The variable in the dataset to filter by, comparing to md_onc_visit_int.  This should be an interval from date of birth in days to make sense.
get_med_onc_by_timing <- function(
  dat_med_onc,
  dat_cutoff,
  var_cutoff,
  remove_missing = T,
  only_keep_latest = T
) {
  if (!("md_ecog_imputed" %in% names(dat_med_onc))) {
    cli_abort(
      "This function assumes you've run augment_med_onc_imputed_ecog(.), run that first"
    )
  }

  if (remove_missing) {
    dat_med_onc %<>% filter(!(md_ecog_imputed %in% "Not documented in note"))
  }

  dat_cutoff %<>% select(record_id, .cut = all_of(var_cutoff))

  rtn <- left_join(
    dat_med_onc,
    dat_cutoff,
    by = "record_id"
  )

  rtn %<>%
    filter(md_onc_visit_int <= (.cut + 0.5)) %>% # tolerance of 0.5 days.
    select(-.cut)

  if (only_keep_latest) {
    rtn %<>%
      group_by(record_id) %>%
      arrange(md_onc_visit_int) %>%
      slice(n()) %>%
      ungroup(.)
  }

  return(rtn)
}
