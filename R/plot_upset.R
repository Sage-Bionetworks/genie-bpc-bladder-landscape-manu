plot_upset <- function(
  dat,
  order_by = "freq",
  y_lab = "n with combo",
  n_set = 15,
  n_intersect = 50
) {
  dat %>%
    as.data.frame %>%
    UpSetR::upset(
      data = .,
      order.by = order_by,
      mainbar.y.label = y_lab,
      nsets = n_set,
      nintersects = n_intersect,
      set_size.show = T
    )
}
