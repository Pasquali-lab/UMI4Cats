#' Create Domainogram
#'
#' Using as input the raw UMIs, this package creates a domainogram for the supplied scales.
#' @inheritParams se
#' @export
#' @return A matrix where the first column represents the fragment end coordinates (start) and the rest represent the
#' number of UMIs found when using a specific scale.
calculateDomainogram <- function(umi4c,
                                 scales=5:150) {

  ## Get cumsum separately upstream and downstream
  cumsum_up <- apply(assay(umi4c)[rowRanges(umi4c)$position=="upstream",], 2, cumsum)
  cumsum_down <- apply(assay(umi4c)[rowRanges(umi4c)$position=="downstream",], 2, cumsum)

  ## Create domainogram
  for(i in 1:ncol(assay(umi4c))) {
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
    rownames(dgram) <- rownames(assay(umi4c))
    colnames(dgram) <- scales

    dgram(umi4c)[[colnames(assay(umi4c))[i]]] <- dgram
  }

  return(umi4c)
}
