#' Test differential contacts
#'
#' Performs fisher's exact tests for each region supplied comparing the umi4c_obj to the reference umi4c_obj,
#' thus testing for differential contacts in such regions.
#' @param umi4c_obj Object as returned by processUMI4C() to use as "treatment".
#' @param ref_umi4c_obj Object as returned by processUMI4C() to use as "reference".
#' @param window_total Region around the bait to use for getting total number of UMIs (the region around the bait
#' is already removed when running processUMI4()). Default: 2Mb.
#' @param padj_method Method to use for adjusting the p-values for multiple testing. See p.adjust().
#' @export
# regions <- raw.gr[floor(runif(n=10, min=1, max=length(raw.gr))),]
# regions <- GenomicRanges::resize(regions, fix="center", width=1e5)
# regions <- regioneR::joinRegions(regions)
differentialContacts <- function(umi4c_obj,
                                 ref_umi4c_obj,
                                 regions,
                                 window_total=1e6,
                                 padj_method="fdr") {
  ## Convert raw reads to GRanges
  raw.gr <- GenomicRanges::GRanges(paste0(GenomicRanges::seqnames(umi4c_obj$bait), ":",
                           umi4c_obj$raw$coord))
  raw.gr$treat <- umi4c_obj$raw$UMIs
  raw.gr$ref <- ref_umi4c_obj$raw$UMIs

  ## Get total in XMb window
  window <- GenomicRanges::resize(umi4c_obj$bait,
                                  fix="center",
                                  width=window_total)

  ## Subset raw by elements in window
  raw.gr <- IRanges::subsetByOverlaps(raw.gr, window)

  ## Get total in window
  total_treat <- sum(raw.gr$treat)
  total_ref <- sum(raw.gr$ref)

  ## Get overlaps with regions of interest
  if (length(unique(GenomicRanges::mcols(regions)))!=length(regions)) {
    regions$ID <- 1:length(regions)
  } else {
    colnames(GenomicRegions::mcols(regions))[1] <- "ID"
  }

  ## Get overlaps of fragments with regions of interest
  ols <- GenomicRanges::findOverlaps(raw.gr, regions)

  raw.ols <- cbind(data.frame(GenomicRanges::mcols(raw.gr))[S4Vectors::queryHits(ols),],
                   data.frame(GenomicRanges::mcols(regions))[S4Vectors::subjectHits(ols),])
  colnames(raw.ols)[3] <- "ID"

  ## Sum UMIs overlapping with each region
  raw.sp <- split(raw.ols, raw.ols$ID)

  raw.sum <- do.call(rbind,
                     lapply(raw.sp,
                            function(x) data.frame(treat=sum(x$treat),
                                                   ref=sum(x$ref),
                                                   ID=unique(x$ID))))

  ## Perform fisher's exact test for each region
  raw.test <- do.call(rbind,
                      lapply(1:nrow(raw.sum),
                             function(x) diffTest(raw.sum$treat[x],
                                                  raw.sum$ref[x],
                                                  total_treat,
                                                  total_ref)))

  ## Adjust pvalue by the number of tests performed
  raw.test$padj <- p.adjust(raw.test$pval, method=padj_method)

  ## Return final data.frame
  final.test <- cbind(raw.sum,
                      raw.test)

  return(final.test)
}


diffTest <- function(treat,
                      ref,
                      total_treat,
                      total_ref) {
  ## Create contingency table
  mat <- matrix(c(treat, ref,
                  total_treat-treat, total_ref-ref),
                ncol=2, byrow=TRUE, dimnames=list(c("inRegion", "notInRegion"),
                                                  c("ctrl", "treat")))

  ## Perform test
  test <- fisher.test(mat)

  ## Export test data.frame
  test.df <- data.frame(pval=test$p.value,
                        OR=test$estimate,
                        CI_lower=test$conf.int[1],
                        CI_upper=test$conf.int[2])

  return(test.df)
}
