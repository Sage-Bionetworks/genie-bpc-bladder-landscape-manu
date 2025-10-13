# Combined summary of mets and stage at dx.
format_stage_mets_dx <- function(
  var_stage_dx,
  var_ca_dmets_yn,
  drop_unused = F
) {
  sm_lev <- c(
    "Stage I",
    "Stage II/III",
    "Stage IV (no met)",
    "Stage IV (met at dx)",
    "Stage IV (met unk.)"
  )

  f <- case_when(
    var_stage_dx %in% c("Stage 0", "Stage I") ~ sm_lev[1],
    # Stage I-III NOS was a clinician decision, obviously could be false.
    var_stage_dx %in% c("Stage II", "Stage III", "Stage I-III NOS") ~ sm_lev[2],
    var_stage_dx %in%
      "Stage IV" &
      var_ca_dmets_yn %in%
        "No - patient is stage IV with no distant metastases" ~
      sm_lev[3],
    var_stage_dx %in% "Stage IV" & var_ca_dmets_yn %in% "Yes" ~ sm_lev[4],
    var_stage_dx %in% "Stage IV" ~ sm_lev[5],
    T ~ NA_character_
  )

  f <- factor(f, levels = sm_lev)

  if (drop_unused) {
    f <- forcats::fct_drop(f)
  }

  return(f)
}
