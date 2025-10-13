sun_wrap <- function(sun_dat, seed) {
  set.seed(seed)
  pal <- make_sun_pal(nrow(sun_dat))
  js <- sunburstR::sunburst(
    sun_dat,
    legend = F,
    count = T,
    colors = pal
  )
  return(js)
}
