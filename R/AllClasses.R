#' UMI4C-class
#'
#' @name UMI4C
#' @aliases UMI4C-class
#' @docType class
#' @note The \code{UMI4C} object extends the \code{DESeqDataSet} class.
#'
#' @param colData Data.frame object containing all the sample information required for the analysis:
#' \enumerate{
#' \item sampleID: Name for the sample.
#' \item replicate: Replicate number.
#' \item condition: Experimental condition or tissue.
#' \item file: Path for the output file containing UMI contact information.
#' }
#'
#' @param metadata List containing the following elements:
#' \enumerate{
#' \item bait: GRanges object representing the position of the bait used for the analysis.
#' }
#'
#' @param assayData Matrix containing the ID for the restriction fragment as row.names and each column representing
#' the number of raw UMIs for a specific sample.
#'
#' @keywords package
#' @rdname UMI4C
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @export
.UMI4C <- setClass("UMI4C",
                   contains = "SummarizedExperiment")

setValidity( "UMI4C", function( object ) {
if (! ("file" %in% names(colData(object))))
    return( "colData must contain 'file'" )
  TRUE
} )

##' @rdname UMI4C
#' @param bait_exclusion Region around the bait (in bp) to be excluded from the analysis. Default: 3kb.
#' @param bait_upstream Number of bp upstream of the bait to use for the analysis. Default: 500kb.
#' @param bait_downstream Number of bp downstream of the bait to use for the analysis. Default: 500kb.
##' @import GenomicRanges SummarizedExperiment
##' @export
UMI4C <- function(colData,
                  viewpoint_name,
                  bait_exclusion=3e3,
                  bait_upstream=5e5,
                  bait_downstream=5e5){
  if (! ("condition" %in% names(colData)))
    return( "colData must contain 'condition'" )
  if (! ("replicate" %in% names(colData)))
    return( "colData must contain 'replicate'" )
  if (! ("file" %in% names(colData)))
    return( "colData must contain 'file'" )

  ## Load UMI4C matrices
  mats <- lapply(as.character(colData$file),
                 read.delim,
                 header=T,
                 stringsAsFactors=F)
  names(mats) <- colData$sampleID

  ## Obtain bait information
  baits <- lapply(mats,
                  function(x) unique(GRanges(seqnames=unique(x[,1]),
                                             ranges=IRanges(start=x[,2],
                                                            end=x[,2]))))
  baits <- unlist(GRangesList(baits))
  bait <- unique(baits)
  if (length(bait) > 1)
    return("Bait position for the supplied samples differ. Please check again your methods.")

  bait$name <- viewpoint_name
  rownames(bait) <- NULL

  ## Create assay matrix with raw UMIs
  umis <- do.call(rbind, mats)
  umis <- umis[,-c(1:2)]
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
  assay <- as.matrix(umis.d[,-c(1:2,ncol(umis.d))])
  rownames(assay) <- rowRanges$id_contact

  ## Create summarizedExperiment
  tmp <- SummarizedExperiment(colData=colData,
                              rowRanges=rowRanges,
                              metadata=list(bait=bait),
                              assays=SimpleList(raw=assay))

  # umi4c <- .UMI4C(tmp)
  UMI4C <- tmp
  ## TODO: findOverlaps() does not work when class is different from SummarizedExperiment

  ## Remove region around bait
  bait_exp <- GenomicRanges::resize(metadata(UMI4C)$bait, fix="center", width=bait_exclusion)
  to_exclude <- subsetByOverlaps(UMI4C, bait_exp)
  UMI4C <- UMI4C[!(rowRanges(UMI4C)$id_contact %in% rowRanges(to_exclude)$id_contact),]

  ## Remove regions outside region
  region <- regioneR::extendRegions(metadata(UMI4C)$bait,
                                    extend.start=bait_upstream,
                                    extend.end=bait_downstream)
  UMI4C <- subsetByOverlaps(UMI4C, region)

  return(UMI4C)
}
