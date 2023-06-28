augment_med_onc_imputed_ecog <- function(med_onc_dat) {
  
  med_onc_dat <- med_onc_dat %>%
    mutate(
      md_ecog_exists = !is.na(md_ecog) & !(md_ecog %in% "Not documented in note"),
      md_karn_exists = !is.na(md_karnof) & !(md_karnof %in% "Not documented in the note")
    )
  
  chk_ecog_karn_exclusive <- med_onc_dat %>% 
    filter(md_ecog_exists & md_karn_exists) %>%
    nrow %>%
    is_in(.,0)
  if (!chk_ecog_karn_exclusive) {
    cli::cli_abort("Karnofsky and ECOG scores both exist in some rows.")
  }
  
  med_onc_dat %<>%
    mutate(md_karn_as_ecog = translate_karn_to_ecog(md_karnof)) %>%
    mutate(
      md_ecog_imputed = case_when(
        !md_ecog_exists & !md_karn_exists ~ "Not documented in note",
        md_ecog_exists ~ md_ecog,
        md_karn_exists ~ md_karn_as_ecog,
        T ~ "Error in imputing ECOG from Karnofsky."
      ),
      md_ecog_imp_source = case_when(
        !md_ecog_exists & !md_karn_exists ~ NA_character_,
        md_ecog_exists ~ "ECOG scale",
        md_karn_exists ~ "Karnofsky (translated)",
        T ~ "Error in imputing ECOG from Karnofsky."
      )
    ) %>%
    select(-c(md_ecog_exists, md_karn_exists))
  
  return(med_onc_dat)
  
}