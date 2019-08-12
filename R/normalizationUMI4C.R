#' Normalize UMI4C
#'
#' Scales the umi4c_obj to the ref_umi4c_obj, creating all necessary elements to produce the comparative plot.
#' @param UMI4C Object outputed by the function \code{\link{processUMI4C}} representing a single UMI-4C processed profile.
#' @param ref_umi4c Object outputed by the function \code{\link{processUMI4C}} representing a single UMI-4C processed profile to use as reference.
#' \code{umi4c_obj} will be scaled to the \code{ref_umi4c_obj}.
#' @param norm_bins Bins used to normalize the number of UMIs from \code{umi4c_obj} to \code{umi4c_ref}. Default: 10^(3:6).
#' @param sd Standard devition value for the adaptative smoothing. Default: 2 (extracted from Schwartzman et al).
#' @export
#' @return A list object with the following named objects:
normalizeUMI4C <- function(UMI4C,
                           ref_umi4c = "smaller",
                           norm_bins=10^(3:6),
                           sd=2) {
  if (ref_umi4c=="smaller") ref_umi4c <- colnames(assay(UMI4C))[which(colSums(assay(UMI4C)) == min(colSums(assay(UMI4C))))]
  
  metadata(UMI4C)$ref_umi4c <- ref_umi4c
  
  UMI4C <- getNormalizationMatrix(UMI4C,
                                  norm_bins=norm_bins)
  
  UMI4C <- normalizeDomainogram(UMI4C)

  ## Plot comparative trend
  UMI4C <- smoothTrendAdaptativeComp(UMI4C,
                                     sd=sd)

  return(UMI4C)
}

#' Get normalization matrix
#'
#' Will return a vector for normalizing domainogram.
#' @inheritParams processComparisonUMI4C
#' @export
getNormalizationMatrix <- function(UMI4C,
                                   norm_bins=10^(3:6),
                                   post_smooth_win=50,
                                   r_expand=1.2) {
  bins <- sort(norm_bins)
  
  bait_start <- GenomicRanges::start(metadata(UMI4C)$bait)
  coords <- start(rowRanges(UMI4C))
  dist <- abs(coords - bait_start)
  
  norm_mat <- matrix(NA,
                     ncol=ncol(assay(UMI4C)),
                     nrow=nrow(assay(UMI4C)))
  
  for (s in 1:ncol(assay(UMI4C))) {
    dgram <- metadata(UMI4C)$dgram[[s]]
    mols <- assay(UMI4C)[,s]
    
    sum_vec <- rep(NA, length(coords))
    
    f_bin <- dist < bins[1]
    tot_in_bin <- sum(mols[f_bin]) # stat=linear
    sum_vec[f_bin] <- tot_in_bin
    
    for (i in 2:length(bins)) {
      f_bin_exp <- dist < bins[i]*r_expand & dist >= bins[i-1]/r_expand
      f_bin <- dist < bins[i] & dist >= bins[i-1]
      
      tot_in_bin <- sum(mols[f_bin_exp])
      sum_vec[f_bin] <- tot_in_bin
      
      message("bin ", bins[i], " total ", tot_in_bin)
    }
    
    f_bin <- dist >= bins[length(bins)]
    tot_in_bin <- sum(mols[f_bin])
    sum_vec[f_bin] <- tot_in_bin
    
    sum_vec <- zoo::rollmean(sum_vec,
                             post_smooth_win,
                             fill=c(sum_vec[1], NA, sum_vec[length(sum_vec)]))
    
    norm_mat[,s] <- sum_vec
  }
  
  colnames(norm_mat) <- colnames(assay(UMI4C))
  rownames(norm_mat) <- rownames(assay(UMI4C))
  
  ## Normalize to reference profile
  norm_vector <- norm_mat[,metadata(UMI4C)$ref_umi4c]
  norm_mat <- norm_vector*(1/norm_mat)
  
  assays(UMI4C)$norm_mat <- norm_mat
  
  return(UMI4C)
}

#' Normalize domainogram to reference
#'
#' Will normalize the domainogram of a given umi4c_obj to the reference.
#' @inheritParams processComparisonUMI4C
#' @export
normalizeDomainogram <- function(UMI4C) {
  
  metadata(UMI4C)$norm_dgram <- list()
  for (s in 1:ncol(assay(UMI4C))) {
    metadata(UMI4C)$norm_dgram[[colnames(assay(UMI4C))[s]]] <- metadata(UMI4C)$dgram[[s]] * assays(UMI4C)$norm_mat[,s]
  }
  
  return(UMI4C)
}


#' Adaptative smoothing of scaled trend
#'
#' Will perform adaptative smoothing will scaling one profile to the reference UMI-4C profile.
#' @inheritParams processComparisonUMI4C
#' @export
smoothTrendAdaptativeComp <- function(UMI4C,
                                      sd=2) {
  bait_coordinates <- metadata(UMI4C)$bait
  
  ## Select min_win_cov
  min_win_cov <- metadata(UMI4C)$min_win_cov[metadata(UMI4C)$ref_umi4c]

  # Get domainograms and convert NA to 0
  dgrams <- metadata(UMI4C)$dgram
  dgrams <- lapply(dgrams,
                   function(x) {x[is.na(x)] <- 0; x})

  # Initialize values
  vector_trend <- matrix(NA, 
                         nrow=nrow(dgrams[[1]]),
                         ncol=length(dgrams))
  scale <- rep(NA, nrow(dgrams[[1]]))
  coord <- rep(NA, nrow(dgrams[[1]]))
  
  base_coord <- start(rowRanges(UMI4C))
  mean_coord <- numeric()
  base_scales <- as.numeric(colnames(dgrams[[1]]))
  
  for (i in 1:ncol(dgrams[[1]])) {
    cur_scale <- as.numeric(colnames(dgrams[[1]])[i])
    s <- rowSums(do.call(cbind, 
                         lapply(dgrams,
                                function(x) x[,i]*cur_scale*2 > min_win_cov)))
    f <- rowSums(is.na(vector_trend))==4 & (s == 4)
    
    ###
    vector_trend[f,] <- do.call(cbind,
                                lapply(dgrams,
                                       function(x) x[f,i]))
    
    ## Get geometric mean coordinates
    mean_coords <- geoMeanCoordinates(coords=base_coord,
                                      scale=cur_scale*2,
                                      bait_start=GenomicRanges::start(metadata(UMI4C)$bait))
    coord[f] <- mean_coords[f]
    scale[f] <- cur_scale
  }
  
  vector_norm <- vector_trend * assays(UMI4C)$norm_mat
  coord[is.na(coord)] <- base_coord[is.na(coord)]
  
  ## Append all info to data.frame
  vector_norm <- data.frame(vector_norm)
  vector_norm$coord <- coord
  vector_norm$scale <- scale
  vector_norm$group <- "upstream"
  vector_norm$group[vector_norm$coord > GenomicRanges::start(metadata(UMI4C)$bait)] <- "downstream"
  
  ## Create long format data.frame
  trend_norm <- reshape2::melt(vector_norm,
                               id.vars=c("coord", "scale", "group"),
                               value.vars=1:(ncol(vector_norm)-3),
                               value.name="trend",
                               variable.name="sample")

  ## Add SD
  dev <- sd * sqrt(trend_norm$trend/(trend_norm$scale*2))
  trend_norm$devP <- trend_norm$trend + dev
  trend_norm$devM <- trend_norm$trend - dev
  
  ## Return comparison object
  metadata(UMI4C)$norm_trend <- trend_norm
  
  return(UMI4C)
}
