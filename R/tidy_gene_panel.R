# Reads in a .txt gene panel file, and returns a tibble with the key info.
tidy_gene_panel <- function(path) {
  raw <- read_yaml(path, readLines.warn = F)
  vec_genes <- stringr::str_split_1(string = raw$gene_list, pattern = "\\t")
  rtn <- tibble(
    stable_id = raw$stable_id,
    hugo = vec_genes,
    tested = T
  )

  return(rtn)
}
