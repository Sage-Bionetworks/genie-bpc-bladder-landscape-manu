
get_lot_most_common_2 <- function(
    dat_lot,
    max_line = 3,
    max_most_common = 5
) {
  
  dat_lot %<>% 
    filter(!is.na(line_therapy)) %>%
    filter(line_therapy <= max_line)
  
  dat_lot %<>%
    group_by(line_therapy) %>%
    mutate(n_with_reg_number = n()) %>%
    ungroup()
  
  dat_lot %<>%
    group_by(line_therapy, regimen_drugs) %>%
    summarize(
      n_drug = n(),
      n_with_reg_number = first(n_with_reg_number),
      .groups = "drop"
    )
  
  dat_lot %<>%
    group_by(line_therapy) %>%
    arrange(desc(n_drug)) %>%
    slice(1:max_most_common) %>%
    mutate(common = paste(as.character(english::ordinal(1:n())),
                          "most common")) %>%
    ungroup
  
  dat_lot %<>%
    arrange(line_therapy, desc(n_drug))
  
  rtn <- dat_lot %>%
    mutate(line_therapy_txt = english::ordinal(line_therapy)) %>%
    select(line_therapy_txt, everything())
  
  return(rtn)
}
