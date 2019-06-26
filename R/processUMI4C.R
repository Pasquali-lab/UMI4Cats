#' Process a single UMI4C profile
#'
#' This function processes a single UMI4C profile and generates all data necessary for plotting and performing differential analysis.
#' @param input A two-column data.frame where the first column represents the fragend and the second the number of UMIs.
#' @param bait_coordinates GRanges object with the coordinates of the bait used for the experiment.
#' @param scales Scales to use for smoothing and plotting the domainogram. Scales indicate how many fragment
#' ends are merged together to obtain the normalized number of UMIs. Default: 5:150.
#' @param bait_exclusion Region around the bait (in bp) to be excluded from the analysis. Default: 3kb.
#' @param bait_upstream Number of bp upstream of the bait to use for the analysis. Default: 500kb.
#' @param bait_downstream Number of bp downstream of the bait to use for the analysis. Default: 500kb.
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
processUMI4C <- function(input,
                         bait_coordinates,
                         scales=5:150,
                         bait_exclusion=3e3,
                         bait_upstream=5e5,
                         bait_downstream=5e5,
                         min_win_cov=NULL,
                         factor=0.025,
                         sd=2) {

  # TODO: Match this section with output of Marc's pipeline

  if (is.character(input)) input <- read.delim(input, stringsAsFactors=FALSE, header=T)

  if (ncol(input)==5) {
    bait_coordinates <- regioneR::toGRanges(unique(input[,1:2]))
    input <- input[,c(4:5)]
  }

  ## Order input by coordinates
  input <- input[order(input[,1]),]

  ## Create domainogram
  dgram <- createDomainogram(input=input,
                             bait_coordinates=bait_coordinates,
                             scales=scales,
                             bait_exclusion=bait_exclusion,
                             bait_upstream=bait_upstream,
                             bait_downstream=bait_downstream)

  ## Get default min_win_cov
  if (is.null(min_win_cov)) min_win_cov <- getMinWinMols(sum(input[,2]),
                                                         factor=factor)

  ## Create adaptative smoothed trend
  trend <- smoothTrendAdaptative(dgram,
                                 min_win_cov=min_win_cov,
                                 bait_coordinates=bait_coordinates,
                                 sd=sd)

  ## Obtain raw
  raw <- filterRaw(input=input,
                   bait_coordinates=bait_coordinates,
                   bait_exclusion=bait_exclusion,
                   bait_upstream=bait_upstream,
                   bait_downstream=bait_downstream)

  ## Create UMI4CObject
  umi4c_obj <- list(raw=raw,
                    bait=bait_coordinates,
                    dgram=dgram,
                    trend=trend)

  return(umi4c_obj)
}

#' Create Domainogram
#'
#' Using as input the raw UMIs, this package creates a domainogram for the supplied scales.
#' @inheritParams processUMI4C
#' @export
#' @return A matrix where the first column represents the fragment end coordinates (start) and the rest represent the
#' number of UMIs found when using a specific scale.
createDomainogram <- function(input,
                              bait_coordinates,
                              scales=5:150,
                              bait_exclusion=3e3,
                              bait_upstream=5e5,
                              bait_downstream=5e5) {

  ## Remove region around bait
  bait_exp <- GenomicRanges::resize(bait_coordinates, fix="center", width=bait_exclusion)
  l <- input[,1] <= GenomicRanges::start(bait_exp) |
    input[,1] >= GenomicRanges::end(bait_exp)
  input <- input[l,]

  ## Remove regions outside region
  region <- regioneR::extendRegions(bait_coordinates,
                                    extend.start=bait_upstream,
                                    extend.end=bait_downstream)

  l <- input[,1] >= GenomicRanges::start(region) &
    input[,1] <= GenomicRanges::end(region)
  input <- input[l,]

  ## Divide upstream and downstream coordintes
  upstream <- input[,1] < GenomicRanges::start(bait_coordinates)
  upstream <- seq(1, nrow(input))[upstream]
  downstream <- input[,1] > GenomicRanges::start(bait_coordinates)
  downstream <- seq(1, nrow(input))[downstream]

  ## Get cumsum separately upstream and downstream
  cumsum_up <- cumsum(input[upstream,2])
  cumsum_down <- cumsum(input[downstream,2])

  ## Create domainogram
  dgram <- matrix(input[,1])
  for (s in scales) {
    dgram_up <- c(rep(NA, s),
                  (cumsum_up[(2*s+1):upstream[length(upstream)]] -
                     cumsum_up[1:(upstream[length(upstream)-2*s])])/(2*s),
                  rep(NA,s))

    dgram_dw <- c(rep(NA, s),
                  (cumsum_down[(2*s+1):length(downstream)] -
                     cumsum_down[1:(length(downstream)-2*s)])/(2*s),
                  rep(NA,s))

    dgram <- cbind(dgram, c(dgram_up, dgram_dw))
  }

  colnames(dgram) <- c("start", scales)

  return(dgram)
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
smoothTrendAdaptative <- function(dgram,
                                  min_win_cov=NULL,
                                  bait_coordinates,
                                  factor=0.025,
                                  sd=2) {
  ## Get default min_win_cov
  if (is.null(min_win_cov)) min_win_cov <- getMinWinMols(sum(input[,2]),
                                                         factor=factor)
  ## Set NAs in dgram to 0
  dgram[is.na(dgram)] <- 0

  ## Define needed parameters
  vector <- rep(NA, nrow(dgram))
  scale <- rep(NA, nrow(dgram))
  coord <- rep(NA, nrow(dgram))

  base_coord <- dgram[,1]
  mean_coord <- numeric()
  base_scales <- as.numeric(colnames(dgram)[-1])

  ## Obtain trend using scales for which number of UMIs is >= min_win_cov
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

  ## Create trend data.frame
  coord[is.na(coord)] <- base_coord[is.na(coord)]
  trend <- data.frame(coord=coord,
                      trend=vector,
                      scale=scale)

  ## Group according if fragment ends are upstream or downstream of the bait
  trend$group <- "upstream"
  trend$group[trend$coord > GenomicRanges::start(bait_coordinates)] <- "downstream"
  trend <- unique(trend)

  ## Add SD
  dev <- sd * sqrt(trend$trend/(trend$scale*2))
  trend$devP <- trend$trend + dev
  trend$devM <- trend$trend - dev

  return(trend)
}

#' Filter raw reads
#'
#' Filters out all fragment ends around the bait and outside the window of interest.
#' @inheritParams processUMI4C
#' @export
#' @return A data.frame with two columns corresponding to the fragment end (\emph{coord}) and the number of
#' UMIs supporting the article at that fragment end (\emph{UMIs}).
filterRaw <- function(input,
                      bait_coordinates,
                      bait_exclusion=3e3,
                      bait_upstream=5e5,
                      bait_downstream=5e5) {
  ## Order input by coordinates
  input <- input[order(input[,1]),]

  ## Remove region around bait
  bait_exp <- GenomicRanges::resize(bait_coordinates, fix="center", width=bait_exclusion)
  l <- input[,1] <= GenomicRanges::start(bait_exp) |
    input[,1] >= GenomicRanges::end(bait_exp)
  input <- input[l,]

  ## Remove regions outside region
  region <- regioneR::extendRegions(bait_coordinates,
                                    extend.start=bait_upstream,
                                    extend.end=bait_downstream)

  l <- input[,1] >= GenomicRanges::start(region) &
    input[,1] <= GenomicRanges::end(region)
  input <- input[l,]

  ## Create raw output
  raw <- as.data.frame(input[,1:2])
  colnames(raw) <- c("coord", "UMIs")

  return(raw)
}

