get_lot_most_common <- function(
  dat_reg,
  max_reg_number = 3,
  max_most_common = 5
) {
  dat_reg <- lot_filter_prep(dat_reg, max_reg = max_reg_number)

  dat_reg %<>%
    group_by(.conseq_reg_number) %>%
    mutate(n_with_reg_number = n()) %>%
    ungroup()

  dat_reg %<>%
    group_by(.conseq_reg_number, regimen_drugs) %>%
    summarize(
      n_drug = n(),
      n_with_reg_number = first(n_with_reg_number),
      .groups = "drop"
    )

  dat_reg %<>%
    group_by(.conseq_reg_number) %>%
    arrange(desc(n_drug)) %>%
    slice(1:max_most_common) %>%
    mutate(
      common = paste(as.character(english::ordinal(1:n())), "most common")
    ) %>%
    ungroup

  dat_reg %<>%
    arrange(.conseq_reg_number, desc(n_drug))

  rtn <- dat_reg %>%
    mutate(line_of_therapy = english::ordinal(.conseq_reg_number)) %>%
    select(-.conseq_reg_number) %>%
    select(line_of_therapy, everything())

  return(rtn)
}
