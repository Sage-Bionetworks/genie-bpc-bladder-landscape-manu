#' @param dat_binary A binary gene matrix created by make_binary_gene_matrix()
gene_fisher_wrapper <- function(
  dat_binary,
  axis_order = NULL
) {
  assoc <- test_fisher_co_occur(
    dat_binary,
    ignore_cols = c("sample_id"),
    top = Inf, # we assume selection is handled already.
    alpha = 0.05,
    axis_order = axis_order
  )

  assoc %<>%
    mutate(
      assoc_lab = glue("{ct_11} / {ct_10+ct_01}"),
      p_value_adj = p.adjust(p.value, method = "BH"),
      cont_table_lab = cont_table_lab_help(
        ct_11 = ct_11,
        ct_01 = ct_01,
        ct_10 = ct_10,
        ct_00 = ct_00
      )
    )

  return(assoc)
}
