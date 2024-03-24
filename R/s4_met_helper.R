#' @title Helper for metastatic classification w.r.t Stage IV patients.
#' 
#' @param dat_vec_stage_dx_iv Vector from the ca_ind dataset.
#' @param dat_vec_ca_dmets_ym Vector from the ca_dmets_yn dataset.
s4_met_helper <- function(
    dat_vec_stage_dx_iv,
    dat_vec_ca_dmets_yn,
    rtn_factor = T
) {
  
  lev_s4_met <- c(
    "Stage IV, with dmets",
    "Stage IV, no dmet noted",
    "(not stage IV)"
  )
  
  rtn_vec <- case_when(
    !(dat_vec_stage_dx_iv %in% "Stage IV") ~ lev_s4_met[3],
    dat_vec_ca_dmets_yn %in% "Yes" ~ lev_s4_met[1],
    T ~ lev_s4_met[2]
  )
  
  rtn_vec <- factor(rtn_vec, levels = lev_s4_met)
  
  return(rtn_vec)
}
