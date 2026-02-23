dt_wrapper <- function(x, font.size = '10pt') {
  DT::datatable(
    x,
    style = "bootstrap4",
    fillContainer = F,
    options = list(
      initComplete = htmlwidgets::JS(
        "function(settings, json) {",
        paste0(
          "$(this.api().table().container()).css({'font-size': '",
          font.size,
          "'});"
        ),
        "}"
      )
    )
  )
}
