
display_lot_most_common <- function(
    dat_lot
) {
  dat_lot %>%
    mutate(
      str = glue("[{round(n_drug/n_with_reg_number*100)}%] {regimen_drugs}")
    ) %>%
    select(
      common, line_therapy_txt, str
    ) %>%
    mutate(
      line_therapy_txt = paste(
        str_to_sentence(line_therapy_txt),
        "LoT"
      )
    ) %>%
    pivot_wider(
      names_from = line_therapy_txt,
      values_from = str
    ) %>%
    mutate(common = str_to_sentence(common))
  
}
