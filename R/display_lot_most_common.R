
display_lot_most_common <- function(
    dat_lot
) {
  dat_lot %>%
    mutate(
      str = glue("[{round(n_drug/n_with_reg_number*100)}%] {regimen_drugs}")
    ) %>%
    select(
      common, line_of_therapy, str
    ) %>%
    mutate(
      line_of_therapy = paste(
        str_to_sentence(line_of_therapy),
        "LoT"
      )
    ) %>%
    pivot_wider(
      names_from = line_of_therapy,
      values_from = str
    ) %>%
    mutate(common = str_to_sentence(common))
  
}