# Returns character vector entries matching a regex pattern
# I use this for finding variables in a dataset (names vary in BPC)
str_filter <- function(vec, pattern, negate = F) {
  str_ind <- stringr::str_which(vec, pattern, negate = negate)
  vec[str_ind]
}
