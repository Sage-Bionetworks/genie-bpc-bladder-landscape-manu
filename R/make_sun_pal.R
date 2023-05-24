make_sun_pal <- function(size) {
  sample(
    c(
      viridisLite::turbo(n = size/2,
                         begin = 0.05, end = 0.35),
      viridisLite::turbo(n = size/2,
                         begin = 0.6, end = 0.95)
    )
  )
}
