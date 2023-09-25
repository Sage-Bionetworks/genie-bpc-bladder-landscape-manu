
format_mutation_effect <- function(vec) {
  lev_mutation_effect <- c(
    "Loss-of-function",
    "Likely Loss-of-function",
    "Unknown",
    "Inconclusive",
    "Likely Neutral",
    "Likely Switch-of-function",
    "Switch-of-function",
    "Likely Gain-of-function",
    "Gain-of-function"
  )
  
  chk <- unique(vec[!(vec %in% lev_mutation_effect)])
  
  if (length(chk) > 0) {
    cli::cli_alert_danger(glue("Uncoded levels exist: {chk}"))
  }
  
  fac <- factor(vec, levels = lev_mutation_effect)
  
  return(fac)
  
}