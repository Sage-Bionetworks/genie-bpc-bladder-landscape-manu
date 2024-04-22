
cont_tab_variable_pair <- function(
    dat,
    var1,
    var2
) {
  dat %>%
    summarize(
      ct_11 = sum(.data[[var1]] %in% 1 & .data[[var2]] %in% 1),
      ct_10 = sum(.data[[var1]] %in% 1 & .data[[var2]] %in% 0),
      ct_01 = sum(.data[[var1]] %in% 0 & .data[[var2]] %in% 1),
      ct_00 = sum(.data[[var1]] %in% 0 & .data[[var2]] %in% 0)
    )
}