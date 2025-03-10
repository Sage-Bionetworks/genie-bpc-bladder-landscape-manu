day_to_year_genie <- function(
    x
) {
  x_name <- deparse(substitute(x))
  # day variables in GENIE empirically have the following two characteristics:
  m <- mean(x, na.rm = T) # usually around 500+, will check conservatively
  s <- sd(x, na.rm = T) # usually 800+, will check a bit more carefully
  if (m < 100) {
    cli_alert_warning("Mean of {x_name} is <100, check that it's truly a day-scale var.")
  } 
  if (s < 500) {
    cli_alert_warning("Mean of {x_name} is <500, check that it's truly a day-scale var.")
  }
  x / 365.25
}