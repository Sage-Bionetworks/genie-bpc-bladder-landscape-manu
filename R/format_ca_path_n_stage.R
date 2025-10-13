# Turns the variable into a factor for easier processing.
format_ca_path_n_stage <- function(col, drop_unused = T) {
  f <- factor(
    col,
    # pulled straight from AACR data guide:
    levels = c("NX", "N0", "N1", "N2", "N3", "N4", "Not Applicable", "Unknown")
  )

  if (drop_unused) {
    f <- forcats::fct_drop(f)
  }

  return(f)
}
