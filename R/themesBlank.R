themeXblank <- function(...) {
  theme <- ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                          axis.title.x=ggplot2::element_blank(),
                          axis.ticks.x=ggplot2::element_blank(),
                          axis.line.x=ggplot2::element_blank(),
                          ...)

  return(theme)
}

themeYblank <- function(...) {
  theme <- ggplot2::theme(axis.text.y=ggplot2::element_blank(),
                          axis.title.y=ggplot2::element_blank(),
                          axis.ticks.y=ggplot2::element_blank(),
                          axis.line.y=ggplot2::element_blank(),
                          ...)

  return(theme)
}

scaleXCoordinates <- function(chr,
                              ...) {
  scale <- ggplot2::scale_x_continuous(name=paste0("Coordinates ", chr, " (Mb)"),
                              labels=function(x) round(x/1e6, 2),
                              ...)
  return(scale)
}
