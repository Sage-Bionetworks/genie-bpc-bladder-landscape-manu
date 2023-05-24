# Turns the variable into a factor for easier processing.
format_ca_dx_how <- function(col, drop_unused = T) {
  f <- factor(
    col,
    # pulled straight from AACR data guide:
    levels = c("Pathology",
               "Imaging",
               "Physical Exam",
               "Other")
    
  )
  
  if (drop_unused) {
    f <- forcats::fct_drop(f)
  }
  
  return(f)
}