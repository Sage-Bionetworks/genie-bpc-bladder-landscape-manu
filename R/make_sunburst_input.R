make_sunburst_input <- function(
  dat,
  var,
  order_var,
  max_depth = NULL,
  remove_dupes = F
) {
  # every step below needs sorted data, so,
  dat <- dat %>%
    arrange(record_id, ca_seq, .data[[order_var]])

  if (remove_dupes) {
    dat <- dat %>%
      group_by(record_id, ca_seq) %>%
      # keep the first row, or any row which is not the same as the previous row.
      filter(1:n() == 1 | .data[[var]] != lag(.data[[var]]))
  }

  if (!is.null(max_depth)) {
    dat <- dat %>%
      group_by(record_id, ca_seq) %>%
      slice(1:max_depth) %>%
      ungroup()
  }

  rtn <- dat %>%
    group_by(record_id, ca_seq) %>%
    summarize(
      path = paste(.data[[var]], collapse = "-"),
      .groups = "drop"
    ) %>%
    count(path) %>%
    arrange(path) %>%
    rename(Pop = n)

  return(rtn)
}
