# Turns the variable into a factor for easier processing.
format_ca_hist_adeno_squamous <- function(col, drop_unused = T) {
  f <- factor(
    col,
    # pulled straight from AACR data guide:
    levels = c("Adenocarcinoma",
               "Carcinoma",
               "Squamous cell",
               "Sarcoma",
               "Small cell carcinoma",
               "Other histologies/mixed tumor")
    
  )
  
  if (drop_unused) {
    f <- forcats::fct_drop(f)
  }
  
  return(f)
}
