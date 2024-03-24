#' @title Helper for metastatic classification w.r.t Stage IV patients.
#' 
#' @param dat_vec_stage_dx_iv Vector from the ca_ind dataset.
#' @param dat_vec_ca_dmets_ym Vector from the ca_dmets_yn dataset.
s4_met_helper <- function(
    dat_vec_stage_dx_iv,
    dat_vec_ca_dmets_yn
) {
  rtn_vec <- case_when(
    !(dat_vec_stage_dx_iv %in% "Stage IV") ~ "(not stage IV)",
    !(dat_vec_ca_dmets_yn %in% "Yes") ~ "Stage IV, with dmets",
    T ~ "Stage IV, no dmet noted"
  )
  
  return(rtn_vec)
}