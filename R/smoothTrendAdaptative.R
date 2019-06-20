smoothTrendAdaptative <- function(dgram,
                                  min_win_cov=50,
                                  bait_coordinates,
                                  sd=2) {
  dgram[is.na(dgram)] <- 0

  vector <- rep(NA, nrow(dgram))
  scale <- rep(NA, nrow(dgram))
  coord <- rep(NA, nrow(dgram))

  base_coord <- dgram[,1]
  mean_coord <- numeric()
  base_scales <- as.numeric(colnames(dgram)[-1])

  for (i in 2:ncol(dgram)) {
    cur_scale <- as.numeric(colnames(dgram)[i])
    f <- is.na(vector) &
      dgram[,i] * cur_scale * 2 > min_win_cov
    vector[f] <- dgram[f,i]

    ## Get geometric mean coordinates
    mean_coords <- geoMeanCoordinates(coords=base_coord,
                                      scale=cur_scale*2,
                                      bait_start=GenomicRanges::start(bait_coordinates))
    coord[f] <- mean_coords[f]
    scale[f] <- cur_scale
  }

  coord[is.na(coord)] <- base_coord[is.na(coord)]
  trend <- data.frame(coord=coord,
                      trend=vector,
                      scale=scale)
  trend$group <- "upstream"
  trend$group[trend$coord > GenomicRanges::start(bait_coordinates)] <- "downstream"
  trend <- unique(trend)

  ## Add SD
  dev <- sd * sqrt(trend$trend/(trend$scale*2))
  trend$devP <- trend$trend + dev
  trend$devM <- trend$trend - dev

  return(trend)
}
