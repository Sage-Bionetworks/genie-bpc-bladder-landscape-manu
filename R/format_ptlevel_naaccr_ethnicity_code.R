# Turns the variable into a factor for easier processing.
format_ptlevel_naaccr_ethnicity_code <- function(col, drop_unused = T) {
  f <- factor(
    col,
    # pulled straight from AACR data guide:
    levels = c(
      "Non-Spanish; non-Hispanic",
      "Mexican (includes Chicano)",
      "Puerto Rican",
      "Cuban",
      "South or Central American (except Brazil)",
      "Other specified Spanish/Hispanic origin (includes European; excludes Dominican Republic)",
      "Spanish NOS or Hispanic NOS or Latino NOS",
      "Spanish surname only",
      "Dominican Republic",
      "Unknown whether Spanish or not"
    )
  )

  f %<>%
    forcats::fct_collapse(
      `Hispanic or Spanish` = c(
        "Mexican (includes Chicano)",
        "Puerto Rican",
        "Cuban",
        "South or Central American (except Brazil)",
        "Other specified Spanish/Hispanic origin (includes European; excludes Dominican Republic)",
        "Spanish NOS or Hispanic NOS or Latino NOS",
        "Spanish surname only",
        "Dominican Republic"
      )
    )

  if (drop_unused) {
    f <- forcats::fct_drop(f)
  }

  return(f)
}
