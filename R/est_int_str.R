#' @title Estimate and (confidence) interval string
#' @description Creates a string representation of 3 numbers:  Point estimate,
#'   lower bound for an interval, upper bound for an interval.  The purpose is
#'   displaying statistical results (prediciton intervals, confidence intervals,
#'   etc).  The default number of digits comes from the NEJM guidelines, which
#'   recommend two digits as a good starting point for association measures
#'   \url{https://www.nejm.org/author-center/new-manuscripts}
#' @param est Estimate, a numeric vector
#' @param lower Lower bound of the interval
#' @param upper Upper bound of the interval
#' @param est_digits Number of digits to round to for the estimate (default =
#'   2).
#' @param ci_digits Number of digits to round to for the CI (lower and upper).
#'   Defaults to est_digits.
#' @param na String to return if estimate is NA (optional - NULL for no action).
#' @param plus_prefix logical.  Add a plus sign to positive estimates?
#' @examples
#' est_int_str(c(1,2), c(3,4), c(5.6, 7.8))
#' est_int_str(c(NA,2), c(3,4), c(5.6, 7.8), est_digits = 0, plus_prefix = TRUE, na = "miss")
#' @export
est_int_str <- function(
  est,
  lower,
  upper,
  est_digits = 2,
  ci_digits = est_digits,
  na = "",
  plus_prefix = T
) {
  est_s <- form_f(est, digits = est_digits)
  lcb_s <- form_f(lower, digits = ci_digits)
  ucb_s <- form_f(upper, digits = ci_digits)

  if (plus_prefix) {
    est_s <- dplyr::if_else(
      est >= 0,
      glue::glue("+{est_s}"),
      glue::glue("{est_s}"),
      glue::glue("{est_s}")
    )
  }

  rtn <- glue::glue("{est_s} ({lcb_s}, {ucb_s})")

  if (!is.null(na)) {
    rtn[is.na(est)] <- glue::glue("{na}")
  }
  return(rtn)
}
