#' UMI4C-class
#'
#' @name UMI4C
#' @aliases UMI4C-class
#' @docType class
#' @note The \code{UMI4C} class extends the \code{DESeqDataSet} class.
#' @slot colData Data.frame containing the information for constructing the UMI4C experiment object. Needs
#' to contain the following columns:
#' \itemize{
#'     \item \code{sampleID}: Unique identifier for the sample.
#'     \item \code{condition}: Condition for performing differential analysis. Can be control and treatment, two different cell types, etc.
#'     \item \code{replicate}: Number for identifying replicates.
#'     \item \code{file}: File as outputed by \code{umi4CatsContacts} function.
#' }
#' @slot rowRanges \code{\link{GenomicRanges}} object with the coordinates for the restriction fragment ends, their IDs and
#' other additional annotation.
#' @slot metadata List containing the following elements:
#' \enumerate{
#'     \item \code{bait}: GRanges object representing the position of the bait used for the analysis.
#'     \item \code{scales}: Numeric vector containing the scales uses for calculating the domainogram.
#'     \item \code{min_win_factor}: Factor for calculating the minimum molecules requiered in a window for not merging it with the
#'     next one when calculating the adaptative smoothing trend.
#'     \item \code{region}: \code{GenomicRanges} object with the coordinates of the genomic window used for analyzing UMI4C data.
#'     \item \code{ref_umi4c}: Name of the sample used as reference for normalization.
#' }
#' @slot assays  Matrix containing the ID for the restriction fragment as row.names and each column a specific sample.
#' After running the \code{makeUMI4C} function will contain the following data:
#' \enumerate{
#'     \item \code{umis}: Raw number of UMIs detected by the \code{UMI4Cats} analysis.
#'     \item \code{norm_mat}: Matrix containing the normalization factors for each sample and fragment end.
#'     \item \code{trend}: Adaptative smoothing trend of UMIs.
#'     \item \code{geo_coords}: Geometric coordinates obtained when performing the adaptative smoothing.
#'     \item \code{scale}: Scale selected for the adaptative smoothing.
#'     \item \code{sd}: Stantard deviation for the adaptative smoothing trend.
#' }
#' @slot dgram List containing the domainograms for each sample. A domainogram is matrix where columns are different scales selected
#' for merging UMI counts and rows are the restriction fragments.
#' @slot results List containing the results for the differential analysis.
#' @rdname UMI4C
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @export
.UMI4C <- setClass("UMI4C",
                   slots = representation(
                     dgram="SimpleList",
                     results="SimpleList"
                   ),
                   contains = "RangedSummarizedExperiment")

setValidity( "UMI4C", function( object ) {
  TRUE
} )

#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
UMI4C <- function(dgram=S4Vectors::SimpleList(),
                  results=S4Vectors::SimpleList(),
                  ...) {
  se <- SummarizedExperiment(...)
  .UMI4C(se,
         dgram=dgram,
         results=results)
}

#' Make UMI4C object
#'
#' The \code{UMI4C-class} constructor is the function \code{makeUMI4C}. By using the arguments listed below,
#' performs the necessary steps to process UMI4C data and summarize it in an object of class \code{UMI4C}
#' @rdname UMI4C
#' @param colData Data.frame containing the information for constructing the UMI4C experiment object. Needs
#' to contain the following columns:
#' \itemize{
#'     \item sampleID Unique identifier for the sample.
#'     \item condition Condition for performing differential analysis. Can be control and treatment, two different cell types, etc.
#'     \item replicate Number for identifying replicates.
#'     \item file File as outputed by \code{umi4CatsContacts} function.
#' }
#' @param viewpoint_name Character indicating the name for the used viewpoint.
#' @param normalized Logical indicating whether UMI4C profiles should be normalized to the sample with less UMIs (ref_umi4c). Default: TRUE
#' @param bait_exclusion Region around the bait (in bp) to be excluded from the analysis. Default: 3kb.
#' @param bait_expansion Number of bp upstream and downstream of the bait to use for the analysis (window centered in
#' bait). Default: 1Mb.
#' @param scales Numeric vector containing the scales for calculating the domainogram.
#' @param min_win_factor Proportion of UMIs that need to be found in a specific window for adaptative trend calculation
#' @param sd Stantard deviation for adaptative trend.
#' @import GenomicRanges
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
makeUMI4C <- function(colData,
                      viewpoint_name="Unknown",
                      normalized=TRUE,
                      bait_exclusion=3e3,
                      bait_expansion=1e6,
                      scales=5:150,
                      min_win_factor=0.02,
                      sd=2){
  if (! ("condition" %in% names(colData)))
    stop( "colData must contain 'condition'" )
  if (! ("replicate" %in% names(colData)))
    stop( "colData must contain 'replicate'" )
  if (! ("file" %in% names(colData)))
    stop( "colData must contain 'file'" )

  colData$sampleID <- gsub(".", "_", colData$sampleID, fixed=TRUE)

  ## Load UMI4C matrices
  mats <- lapply(as.character(colData$file),
                 read.delim,
                 header=T,
                 stringsAsFactors=F)
  names(mats) <- colData$sampleID

  nrows <- sapply(mats, nrow)
  if (length(unique(nrows))!=1) stop("Number of rows for the supplied files is different. Please check again your methods.")

  pos <- lapply(mats, function(x) paste0(x[,3], ":", x[,4]))
  if (length(unique(pos))!=1) stop("Fragment end coordinates for your files are different. Please check again your methods.")

  ## Obtain bait information
  baits <- lapply(mats,
                  function(x) unique(GRanges(seqnames=unique(x[,1]),
                                             ranges=IRanges(start=x[,2],
                                                            end=x[,2]))))
  baits <- unlist(GRangesList(baits))
  bait <- unique(baits)
  if (length(bait) > 1)
    stop("Bait position for the supplied samples differ. Please check again your methods.")

  bait$name <- viewpoint_name
  rownames(bait) <- NULL

  ## Create assay matrix with raw UMIs
  umis <- do.call(rbind, mats)
  umis <- umis[,-c(1:2)]
  colnames(umis) <- c("chr_contact", "pos_contact", "UMIs")
  umis$sampleID <- unlist(lapply(strsplit(rownames(umis), ".", fixed=T),
                                 function(x) x[1]))
  umis.d <- reshape2::dcast(umis,
                            chr_contact+pos_contact~sampleID,
                            value.var="UMIs")

  umis.d <- umis.d[order(umis.d$chr_contact, umis.d$pos_contact),]
  umis.d$id_contact <- paste0("frag_", 1:nrow(umis.d))

  # Create row_ranges
  rowRanges <- GRanges(seqnames=umis.d$chr_contact,
                       ranges=IRanges(start=umis.d$pos_contact,
                                      end=umis.d$pos_contact))
  rowRanges$id_contact <- umis.d$id_contact

  # Create assay matrix
  assay <- as.matrix(umis.d[,-c(1:2,ncol(umis.d)),drop=F], labels=TRUE)
  rownames(assay) <- rowRanges$id_contact

  ## Create summarizedExperiment
  umi4c <- UMI4C(colData=colData,
                 rowRanges=rowRanges,
                 metadata=list(bait=bait,
                               scales=scales,
                               min_win_factor=min_win_factor),
                 assays=SimpleList(umis=assay))

  ## Remove region around bait
  bait_exp <- GenomicRanges::resize(metadata(umi4c)$bait,
                                    fix="center",
                                    width=bait_exclusion)
  umi4c <- subsetByOverlaps(umi4c, bait_exp, invert=TRUE)

  ## Remove regions outside scope
  region <- GenomicRanges::resize(metadata(umi4c)$bait,
                                  fix="center",
                                  width=bait_expansion)
  umi4c <- subsetByOverlaps(umi4c, region)
  metadata(umi4c)$region <- region

  ## Divide upstream and downstream coordintes
  rowRanges(umi4c)$position <- NA
  rowRanges(umi4c)$position[start(rowRanges(umi4c)) < start(metadata(umi4c)$bait)] <- "upstream"
  rowRanges(umi4c)$position[start(rowRanges(umi4c)) > start(metadata(umi4c)$bait)] <- "downstream"

  ## Get normalization matrix
  metadata(umi4c)$ref_umi4c <- colnames(assay(umi4c))[which(colSums(assay(umi4c))==min(colSums(assay(umi4c))))]
  umi4c <- getNormalizationMatrix(umi4c)

  ## Calculate domainograms
  umi4c <- calculateDomainogram(umi4c,
                                scales=scales,
                                normalized=normalized)

  ## Calculate adaptative trend
  umi4c <- calculateAdaptativeTrend(umi4c,
                                    sd=sd,
                                    normalized=normalized)

  return(umi4c)
}
