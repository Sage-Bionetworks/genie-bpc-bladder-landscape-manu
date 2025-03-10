
# Lines of therapy categories, originally for the 2025 ASCO poster.
categorize_lines <- function(reg_vec, factorize = T) {
  reg_vec <- tolower(reg_vec)
  
  lot_cat_lev <- c(
    'investigational (masked)',
    'immunotherapy',
    'cisplatin-based',
    'carboplatin-based',
    'taxane monotherapy',
    'pemetrexed monotherapy',
    'other'
  )
  
  rtn <- case_when(
    # after a long discussion we decided to include nivolumab here even though
    #   it's not in the 2L io figure.
    str_detect(reg_vec, 'investigation') ~ lot_cat_lev[1],
    str_detect(reg_vec, 'atezoliz|pembroliz|nivolum') ~ lot_cat_lev[2],
    str_detect(reg_vec, 'cisplatin') ~ lot_cat_lev[3],
    str_detect(reg_vec, 'carboplatin') ~ lot_cat_lev[4],
    reg_vec %in% c('docetaxel', 'paclitaxel') ~ lot_cat_lev[5],
    reg_vec %in% c('pemetrex') ~ lot_cat_lev[6],
    T ~ lot_cat_lev[7]
  )
  
  if (factorize) { 
    rtn <- factor(rtn, levels = lot_cat_lev)
  }
  
  return(rtn)
}
