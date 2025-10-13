# Turns the variable into a factor for easier processing.
format_stage_dx <- function(col, drop_unused = T) {
  f <- factor(
    col,
    # pulled straight from AACR data guide:
    levels = c(
      "Stage 0",
      "Stage I",
      "Stage II",
      "Stage III",
      "Stage I-III NOS",
      "Stage IV"
    )
  )

  if (drop_unused) {
    f <- forcats::fct_drop(f)
  }

  return(f)
}
