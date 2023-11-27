# It is assumed this has already been given a column 'var_class'.
assign_lot <- function(
    dat_reg
) {
  if (!('var_class' %in% colnames(dft_reg_dmet_s))) {
    dat_reg %<>% 
      mutate(var_class = NA_character_)
  }
  
  dat_reg %<>%
    arrange(record_id, ca_seq, regimen_number)
  
  dat_reg %<>%
    group_by(record_id, ca_seq) %>%
    mutate(
      line_therapy = lot_helper_one_case(
        var_class = var_class
      )
    ) %>%
    ungroup(.)
  
  return(dat_reg)
  
}


lot_helper_one_case <- function(var_class) {
  lot <- c()
  observed_var_classes <- c()
  lot_iter <- 0L
  # for loop horror:
  for (i in 1:length(var_class)) {
    if (var_class[i] %in% observed_var_classes) {
      lot[i] <- NA_real_
    } else {
      lot_iter <- lot_iter + 1L
      lot[i] <- lot_iter
      if (!is.na(var_class[i])) { 
        observed_var_classes <- c(observed_var_classes, var_class[i])
      }
    }
  }
  return(lot)
}
