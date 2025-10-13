km_tidy_no_covariate <- function(
  dat,
  event_time,
  event_ind,
  trunc_time = NULL
) {
  if (is.null(trunc)) {
    surv_obj <- Surv(time = dat[[event_time]], event = dat[[event_ind]])
  } else {
    surv_obj <- with(
      dat,
      Surv(
        time = dat[[trunc_time]],
        time2 = dat[[event_time]],
        event = dat[[event_ind]]
      )
    )
  }
  survfit(surv_obj ~ 1, data = dat) %>%
    broom::glance(.)
}
