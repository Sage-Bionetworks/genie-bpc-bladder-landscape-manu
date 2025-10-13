n_pct_str <- function(n, d, show_d = F, digits = 0, na = "") {
  if (any((n %% 1 != 0) | (d %% 1 != 0), na.rm = T)) {
    stop("n, d must be integers (n %% 1 == 0 must be true)")
  }
  if (length(na) > 1) {
    stop("na must be a length 1 vector (or NULL)")
  }

  pct <- n / d
  if (show_d) {
    rtn <- glue::glue("{n}/{d} ({form_f(pct*100, digits = digits)}%)")
  } else {
    rtn <- glue::glue("{n} ({form_f(pct*100, digits = digits)}%)")
  }

  if (!is.null(na)) {
    rtn[is.na(n)] <- glue::glue("{na}")
  }
  return(rtn)
}
