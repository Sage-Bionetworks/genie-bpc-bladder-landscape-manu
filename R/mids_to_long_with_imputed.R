# This is a helper for plot_gg_strip - which makes a ggplot version of a diagnostic plot.
#' @param mids_obj A "multiple imputed data set" object created with mice()
#' @param var A character vector the the name of the variable to visualize.
mids_to_long_with_imputed <- function(
    mids_obj,
    var
) {
  rtn <- mice::complete(mids_obj, action = "long", include = T) %>%
    group_by(.id) %>%
    mutate(.imp = as.integer(.imp),
           .imputed = !is.na(.data[[var]]) & is.na(.data[[var]][.imp %in% 0])) %>%
    ungroup(.)
  
  return(rtn)
}
