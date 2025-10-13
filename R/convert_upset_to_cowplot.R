convert_upset_to_cowplot <- function(up_plot, rel_heights = c(3, 2)) {
  # Note:  If you ever want the size bar back, you can do a 4 panel grid with
  #. a NULL upper left plot.
  cp <- cowplot::plot_grid(
    up_plot$Main_bar,
    up_plot$Matrix,
    nrow = 2,
    align = 'hv',
    rel_heights = rel_heights
  )
}
