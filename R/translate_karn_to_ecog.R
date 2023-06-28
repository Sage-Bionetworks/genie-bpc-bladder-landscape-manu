
# Translate Karnofsky scores to ECOG scores.
# Source:  https://oncologypro.esmo.org/oncology-in-practice/practice-tools/performance-scales

translate_karn_to_ecog <- function(col) {
  ecog_levs <- c(
    # pulled from the data:
    "0: Fully active, able to carry on all pre-disease performance without restriction",
    "1: Restricted in physically strenuous activity but ambulatory and able to carry out work of a light or sedentary nature, e.g., light house work, office work",
    "2: Ambulatory and capable of all selfcare but unable to carry out any work activities; up and about more than 50% of waking hours",
    "3: Capable of only limited selfcare; confined to bed or chair more than 50% of waking hours",
    "4: Completely disabled; cannot carry on any selfcare; totally confined to bed or chair",
    "Not documented in note"           
  )
  
  # The Karnofsky strings will just be the first three characters.
  translation_tab <- tribble(
    ~karn_str, ~ecog_str,
    "10:", ecog_levs[5], # ECOG grade of 4.
    "20:", ecog_levs[5],
    "30:", ecog_levs[4],
    "40:", ecog_levs[4],
    "50:", ecog_levs[3],
    "60:", ecog_levs[3],
    "70:", ecog_levs[2],
    "80:", ecog_levs[2],
    "90:", ecog_levs[1], # ECOG grade of 0
    "100", ecog_levs[1],
    "Not", ecog_levs[6]
  )
  
  # t for truncated.
  col_t <- stringr::str_sub(col, 1, 3)
  
  if (any(!(col_t %in% translation_tab$karn_str | is.na(col_t)))) {
    cli::cli_alert_danger(text = "Unexpected strings in 'col'.")
  }
  
  trans_col <- tibble(karn_str = col_t) %>%
    left_join(., translation_tab, by = "karn_str") %>%
    pull(ecog_str)
  
  return(trans_col)
  
}
