
#' @param dat_alt The alterations dataset.
#' @param vec_sample The vector of samples we want included.  This allows us to
#'   include samples which have no alterations found.
make_binary_gene_matrix <- function(dat_alt, vec_sample) {
  rtn <- dat_alt %>%
    filter(sample_id %in% vec_sample) %>%
    group_by(sample_id, hugo) %>%
    summarize(exists = n() >= 1, .groups = "drop") %>%
    pivot_wider(
      names_from = hugo,
      values_from = exists
    )
  
  # Add in rows where nothing at all was found:
  rtn <- left_join(
    tibble(sample_id = vec_sample),
    rtn,
    by = "sample_id"
  )
  
  rtn %<>%
    mutate(
      across(
        .cols = -sample_id,
        .fns = (function(x) {
          x <- as.integer(x)
          x <- if_else(is.na(x), 0L, x)
        })
      )
    ) 
  return(rtn)
  
}