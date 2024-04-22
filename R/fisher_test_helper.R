
fisher_test_helper <- function(
    cell_11, cell_10, cell_01, cell_00,
    alpha = 0.05,
    keep_cols = c('estimate', 'p.value', 'conf.low', 'conf.high')
) {
  rtn <- fisher.test(
    x = matrix(
      c(cell_11, cell_10, cell_01, cell_00),
      nrow = 2, byrow = T
    ),
    conf.level = 1-alpha
  )
  
  rtn %<>%
    broom::tidy(.) %>%
    select(all_of(keep_cols))
  
  return(rtn)
}