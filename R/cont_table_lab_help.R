cont_table_lab_help <- function(
  ct_11,
  ct_01,
  ct_10,
  ct_00,
  width = 5
) {
  sp_help <- function(x) {
    str_pad(x, side = "right", pad = " ", width = width)
  }

  ct_11 <- sp_help(ct_11)
  ct_10 <- sp_help(ct_10)
  ct_01 <- sp_help(ct_01)
  ct_00 <- sp_help(ct_00)

  rtn <- glue("{ct_11}{ct_01}\n{ct_10}{ct_00}")

  return(rtn)
}
