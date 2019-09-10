#' Adaptative smoothing of scaled trend
#'
#' Will perform adaptative smoothing will scaling one profile to the reference UMI-4C profile.
#' @param umi4c umi4c object.
#' @param sd Standard deviation
#' @param scaled Logical indicating whether smoothen trends should be scaled to the sample with less UMIs.
#' @export
calculateAdaptativeTrend <- function(umi4c,
                                     sd=2,
                                     scaled=FALSE) {
  ## Select min_win_cov
  mols <- colSums(assay(umi4c))
  min_mols <- mols * metadata(umi4c)$min_win_factor

  if (scaled) min_mols[min_mols!=min(min_mols)] <- min(min_mols)

  ## Get domainograms and convert NA to 0
  dgrams <- dgram(umi4c)
  dgrams <- lapply(dgrams,
                   function(x) {x[is.na(x)] <- 0; x})

  # Initialize values
  vector_trend <- matrix(NA,
                         nrow=nrow(dgrams[[1]]),
                         ncol=length(dgrams))

  if (scaled) {
    scale <- rep(NA, nrow(dgrams[[1]]))
    coord <- rep(NA, nrow(dgrams[[1]]))
  } else {
    scale <- matrix(NA, nrow=nrow(dgrams[[1]]),
                    ncol=length(dgrams))
    coord <- matrix(NA, nrow=nrow(dgrams[[1]]),
                    ncol=length(dgrams))
  }

  base_coord <- start(rowRanges(umi4c))
  mean_coord <- numeric()

  for (i in 1:ncol(dgrams[[1]])) { # Loop for each scale in dgram
    current_scale <- metadata(umi4c)$scales[i]

    # Get geometric mean coordinates
    mean_coords <- geoMeanCoordinates(coords=base_coord,
                                      scale=current_scale*2,
                                      bait_start=GenomicRanges::start(metadata(umi4c)$bait))

    # Obtain rows in which UMIs > threshold for each sample
    pass <- mapply(function(x, min_win_cov) x[,i]*current_scale*2 > min_win_cov,
                   x=dgrams,
                   min_win_cov=min_mols)

    if (scaled) {
      # Select regions that are still NA in trend and pass cov threshold on all samples
            sums <- rowSums(pass)
      sel <- rowSums(is.na(vector_trend))==length(dgrams) & (sums == length(dgrams))

      vector_trend[sel,] <- do.call(cbind,
                                  lapply(dgrams,
                                         function(x) x[sel,i]))
      coord[sel] <- mean_coords[sel]
      scale[sel] <- current_scale
    } else {
      tmp_dgram <- do.call(cbind, lapply(dgrams, function(x) x[,i]))

      # Select regions that are still NA in trend and pass cov threshold, sample-wise
      na <- is.na(vector_trend)
      pass <- (na + pass) == 2

      vector_trend[pass] <- tmp_dgram[pass]

      mean_coords_mat <- matrix(rep(mean_coords, length(dgrams)), ncol=length(dgrams))
      coord[pass] <- mean_coords_mat[pass]
      scale[pass] <- current_scale
    }
  }

  if(!is.null(assays(umi4c)$norm_mat) & scaled)   vector_trend <- vector_trend * assays(umi4c)$norm_mat

  if (scaled) {
    coord <- matrix(rep(coord, length(dgrams)), ncol=length(dgrams))
    scale <- matrix(rep(scale, length(dgrams)), ncol=length(dgrams))
  }

  base_coord_mat <- matrix(rep(base_coord, length(dgrams)), ncol=length(dgrams))
  coord[is.na(coord)] <- base_coord_mat[is.na(coord)]

  assays(umi4c)$trend <- vector_trend
  assays(umi4c)$geo_coords <- coord
  assays(umi4c)$scale <- scale

  ## Add SD
  dev <- sd * sqrt(vector_trend/(scale*2))
  assays(umi4c)$sd <- dev

  return(umi4c)
}
