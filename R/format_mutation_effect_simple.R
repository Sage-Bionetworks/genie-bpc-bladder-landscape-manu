
format_mutation_effect_simple <- function(vec) {
  lev_me <- c(
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
  
  chk <- unique(vec[!(vec %in% lev_me)])
  
  if (length(chk) > 0) {
    cli::cli_alert_danger(glue("Uncoded levels exist: {chk}"))
  }
  
  fac <- case_when(
    vec %in% lev_me[1:2] ~ "Loss",
    vec %in% lev_me[8:9] ~ "Gain",
    vec %in% lev_me[6:7] ~ "Switch",
    T ~ "Other"
  ) 
  
  fac <- factor(fac, levels = c("Loss", "Other", "Switch", "Gain"))
  
  return(fac)
  
}
