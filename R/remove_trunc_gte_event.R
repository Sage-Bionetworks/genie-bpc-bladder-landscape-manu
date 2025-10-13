remove_trunc_gte_event <- function(
  dat,
  trunc_var,
  event_var,
  message_rows_removed = T
) {
  rtn <- dat %>%
    filter(
      .data[[trunc_var]] < .data[[event_var]]
    )

  if (message_rows_removed) {
    nrow_removed <- nrow(dat) - nrow(rtn)
    cli::cli_alert_info(
      "{nrow_removed} rows removed due to {trunc_var} >= {event_var}."
    )
  }
  return(rtn)
}
