# Turns the variable into a factor for easier processing.
format_md_ecog <- function(col, drop_unused = F) {
  f <- factor(
    str_sub(col, end = 20L), # just to trim this mess down.
    # full level explanation in AACR data guide
    levels = c(
      "0: Fully active, abl",
      "1: Restricted in phy",
      "2: Ambulatory and ca",
      "3: Capable of only l",
      "4: Completely disabl"
    ),
    # Made up by me to be reasonably succinct:
    labels = c(
      "0: Fully active",
      "1: Restricted",
      "2: Ambulatory/selfcare",
      "3: Limited selfcare",
      "4: Completely disabled"
    )
  )
  
  if (drop_unused) {
    f <- forcats::fct_drop(f)
  }
  
  return(f)
}
