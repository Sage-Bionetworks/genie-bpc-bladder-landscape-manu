# Turns the variable into a factor for easier processing.
format_ca_dmets_yn <- function(col, drop_unused = F) {
  f_levs <- c(
    "Stage IV, with dmets",
    "Stage IV, no dmet",
    "Stage IV, unknown",
    "Not Stage IV"
  )

  # Slight rephrase on these:
  col <- case_when(
    col %in% "Yes" ~ f_levs[1],
    col %in% "No - patient is stage IV with no distant metastases" ~ f_levs[2],
    col %in% "Unknown or Not mentioned" ~ f_levs[3],
    is.na(col) ~ f_levs[4],
    T ~ "error"
  )

  f <- factor(
    col,
    levels = f_levs
  )

  if (drop_unused) {
    f <- forcats::fct_drop(f)
  }

  return(f)
}
