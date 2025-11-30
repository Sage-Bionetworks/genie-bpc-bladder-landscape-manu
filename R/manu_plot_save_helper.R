manu_plot_save_helper <- function(
  plot,
  dir,
  name,
  width,
  height,
  as_pdf = T,
  as_png = T,
  as_jpeg = T
) {
  fs::dir_create(dir)
  if (as_pdf) {
    ggsave(
      plot,
      height = height,
      width = width,
      filename = here(dir, paste0(name, '.pdf'))
    )
  }
  if (as_png) {
    ggsave(
      plot,
      height = height,
      width = width,
      filename = here(dir, paste0(name, '.png')),
      dpi = 300
    )
  }
  if (as_jpeg) {
    ggsave(
      plot,
      height = height,
      width = width,
      filename = here(dir, paste0(name, '.jpeg')),
      dpi = 300
    )
  }

  invisible()
}
