#' Process a single UMI4C profile
#'
#' This function processes a single UMI4C profile and generates all data necessary for plotting and performing differential analysis.
#' @param UMI4C UMI4C objected created with the \code{\link{UMI4C}} function.
#' @param scales Scales to use for smoothing and plotting the domainogram. Scales indicate how many fragment
#' ends are merged together to obtain the normalized number of UMIs. Default: 5:150.
#' @param min_win_cov Minimum number of UMIs requiered in a region to perform the adaptaive smoothing. Default: NULL,
#' number is calculated by the function \code{\link{getMinWinMols}} using the factor specified below.
#' @param factor Factor representing the proportion of UMIs from the total UMIs in the analyzed region needed as a minimum
#' to perform the adaptative smoothing.
#' @param sd Standard devition value for the adaptative smoothing. Default: 2 (extracted from Schwartzman et al).
#' @export
#' @return This function returns a list object with the following fields:
#' \itemize{
#'  \item \strong{raw}: A data.frame including the raw UMIs in each fragment end considered for the analysis
#'  (see \code{\link{filterRaw}}).
#'  \item \strong{bait}: A GRanges object containing the coordinates of the bait used in the UMI-4C assay.
#'  \item \strong{dgram}: A matrix containing the domainogram (see \code{\link{createDomainogram}}).
#'  \item \strong{trend}: A data.frame containing the adaptative smoothen trend for the assay
#'  (see \code{\link{smoothTrendAdaptative}}).
#'  }
processUMI4C <- function(UMI4C,
                         scales=5:150,
                         min_win_cov=NULL,
                         factor=0.025,
                         sd=2) {

  ## Create domainogram
  UMI4C <- createDomainogram(UMI4C=UMI4C,
                             scales=scales)

  ## Get default min_win_cov
  if (is.null(min_win_cov)) min_win_cov <- sapply(colSums(assay(UMI4C)),
                                                  getMinWinMols,
                                                  factor=factor)
  metadata(UMI4C)$min_win_cov <- min_win_cov

  ## Create adaptative smoothed trend
  UMI4C <- smoothTrendAdaptative(UMI4C,
                                 sd=sd)

  return(UMI4C)
}

#' Create Domainogram
#'
#' Using as input the raw UMIs, this package creates a domainogram for the supplied scales.
#' @inheritParams processUMI4C
#' @export
#' @return A matrix where the first column represents the fragment end coordinates (start) and the rest represent the
#' number of UMIs found when using a specific scale.
createDomainogram <- function(UMI4C,
                              scales=5:150) {
  ## Divide upstream and downstream coordintes
  rowRanges(UMI4C)$position <- NA
  rowRanges(UMI4C)$position[start(rowRanges(UMI4C)) < start(metadata(UMI4C)$bait)] <- "upstream"
  rowRanges(UMI4C)$position[start(rowRanges(UMI4C)) > start(metadata(UMI4C)$bait)] <- "downstream"

  ## Get cumsum separately upstream and downstream
  cumsum_up <- apply(assay(UMI4C)[rowRanges(UMI4C)$position=="upstream",], 2, cumsum)
  cumsum_down <- apply(assay(UMI4C)[rowRanges(UMI4C)$position=="downstream",], 2, cumsum)

  ## TODO: Understand how to integrate dgram in SummarizedExperiment object
  ## Integrate in metadata

  metadata(UMI4C)$dgram <- list()

  ## Create domainogram
  for(i in 1:ncol(assay(UMI4C))) {
    dgram_up <- matrix(NA, nrow=nrow(cumsum_up), ncol=length(scales))
    dgram_dw <- matrix(NA, nrow=nrow(cumsum_down), ncol=length(scales))

    min_scale <- scales[1]
    for (s in scales) {
      dgram_up[,s-(min_scale-1)] <- c(rep(NA, s),
                         (cumsum_up[(2*s+1):nrow(cumsum_up),i] -
                            cumsum_up[1:(nrow(cumsum_up)-2*s),i])/(2*s),
                         rep(NA,s))

       dgram_dw[,s-(min_scale-1)] <- c(rep(NA, s),
                         (cumsum_down[(2*s+1):nrow(cumsum_down),i] -
                            cumsum_down[1:(nrow(cumsum_down)-2*s),i])/(2*s),
                         rep(NA,s))
      }

      dgram <- rbind(dgram_up,
                     dgram_dw)
      rownames(dgram) <- rownames(assay(UMI4C))
      colnames(dgram) <- scales

      metadata(UMI4C)$dgram[[colnames(assay(UMI4C))[i]]] <- dgram
  }

  return(UMI4C)
}

#' Adaptative smoothing of the UMI-4C trend
#'
#' Performs adaptative smoothing to plot a UMI-4C contact trend. This means that joins fragment ends together until the
#' minimum number of UMIs in the joined fragment is >= min_win_cov.
#' @param dgram Output of the function \code{\link{createDomainogram}}.
#' @inheritParams processUMI4C
#' @export
#' @return A data.frame containing the necessary information to draw the different elements of the UMI-4C trend:
#' \itemize{
#' \item \strong{coord}: Coordinates of the fragment end.
#' \item \strong{trend}: Smoothed UMIs value according to min_win_cov.
#' \item \strong{scale}: Scale used (# of fragends to merge) to achieve min_win_cov.
#' \item \strong{group}: Fragment ends grouped according to postion using bait as reference (upstream or downstream).
#' \item \strong{devP}: Upper bound of standard deviation.
#' \item \strong{devM}: Lower bound of standard deviation.
#' }
smoothTrendAdaptative <- function(UMI4C,
                                  sd=2) {

  trend.all <- data.frame()

  for (s in 1:ncol(assay(UMI4C))) {
    dgram <- metadata(UMI4C)$dgram[[s]]

    ## Set NAs in dgram to 0
    dgram[is.na(dgram)] <- 0

    ## Define needed parameters
    vector <- rep(NA, nrow(dgram))
    scale <- rep(NA, nrow(dgram))
    coord <- rep(NA, nrow(dgram))

    base_coord <- GenomicRanges::start(rowRanges(UMI4C))
    mean_coord <- numeric()
    base_scales <- as.numeric(colnames(dgram))

    ## Select min_win_cov
    if (length(min_win_cov>1)) {
      sel_min_win_cov <- min_win_cov[s]
    } else {
      sel_min_win_cov <- min_win_cov
    }

    ## Obtain trend using scales for which number of UMIs is >= min_win_cov
    for (i in 1:ncol(dgram)) {
      cur_scale <- as.numeric(colnames(dgram)[i])
      f <- is.na(vector) &
        dgram[,i] * cur_scale * 2 > sel_min_win_cov
      vector[f] <- dgram[f,i]

      ## Get geometric mean coordinates
      mean_coords <- geoMeanCoordinates(coords=base_coord,
                                        scale=cur_scale*2,
                                        bait_start=GenomicRanges::start(metadata(UMI4C)$bait))
      coord[f] <- mean_coords[f]
      scale[f] <- cur_scale
    }

    ## Create trend data.frame
    coord[is.na(coord)] <- base_coord[is.na(coord)]
    trend <- data.frame(coord=coord,
                        trend=vector,
                        scale=scale)

    ## Group according if fragment ends are upstream or downstream of the bait
    trend$group <- "upstream"
    trend$group[trend$coord > GenomicRanges::start(metadata(UMI4C)$bait)] <- "downstream"
    trend <- unique(trend)

    ## Add SD
    dev <- sd * sqrt(trend$trend/(trend$scale*2))
    trend$devP <- trend$trend + dev
    trend$devM <- trend$trend - dev

    ## Add sample name
    trend$sample <- colnames(assay(UMI4C))[s]

    trend.all <- rbind(trend.all, trend)
  }

  metadata(UMI4C)$trend <- trend.all

  return(UMI4C)
}
