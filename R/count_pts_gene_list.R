count_pts_gene_list <- function(
    sample_df,
    alt_df,
    gene_list
) {
  if(any(duplicated(sample_df$patient_id))) {
    cli_abort("Duplicate samples detected - won't work for that")
  }
  
  rtn <- alt_df %>%
    filter(hugo %in% gene_list) %>%
    group_by(sample_id) %>%
    summarize(any_alt = 1)
  
  rtn <- sample_df %>%
    select(patient_id, sample_id) %>%
    left_join(., rtn, by = 'sample_id') %>%
    replace_na(list(any_alt = 0))
  
  return(rtn)
  
}