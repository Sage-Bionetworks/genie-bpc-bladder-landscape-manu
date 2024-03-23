collapse_site_mets <- function(
    dat,
    pattern = "^dx_to_dmets_.+_yrs$",
    rtn_col_name = "dx_to_dmets_yrs"
) {
  n_cols_pattern <- sum(str_detect(colnames(dat), pattern), na.rm = T)
  if (n_cols_pattern < 1) {
    cli::cli_abort("No columns match that pattern")
  }
  
  rtn <- dat %>%
    select(record_id, ca_seq, matches(pattern)) %>%
    pivot_longer(
      cols = -c(record_id, ca_seq)
    )
  
  rtn %<>%
    group_by(record_id, ca_seq) %>%
    summarize(
      {{rtn_col_name}} := suppressWarnings(min(value, na.rm = T)),
      .groups = "drop"
    ) %>%
    mutate(
      {{rtn_col_name}} := if_else(
        is.infinite(.data[[rtn_col_name]]), 
        NA_real_, 
        .data[[rtn_col_name]]
      )
    )
  
  return(rtn)
}
