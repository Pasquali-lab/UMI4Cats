#' UMI4C-class
#'
#' @name UMI4C
#' @aliases UMI4C-class
#' @docType class
#' @note The \code{UMI4C} class extends the \linkS4class{SummarizedExperiment} class.
#' @slot colData Data.frame containing the information for constructing the
#' UMI4C experiment object. Needs to contain the following columns:
#' \itemize{
#'     \item \code{sampleID}: Unique identifier for the sample.
#'     \item \code{condition}: Condition for performing differential analysis.
#'     Can be control and treatment, two different cell types, etc.
#'     \item \code{replicate}: Number or ID for identifying different replicates.
#'     \item \code{file}: Path to the files outputed by \code{\link{contactsUMI4C}}.
#' }
#' @slot rowRanges \linkS4class{GRanges} object with the coordinates for
#' the restriction fragment ends, their IDs and other additional annotation columns.
#' @slot metadata List containing the following elements:
#' \enumerate{
#'     \item \code{bait}: \linkS4class{GRanges} object representing the position
#'      of the bait used for the analysis.
#'     \item \code{scales}: Numeric vector containing the scales used for
#'     calculating the domainogram.
#'     \item \code{min_win_factor}: Factor for calculating the minimum molecules
#'      requiered in a window for not merging it with the next one when
#'      calculating the adaptative smoothing trend.
#'     \item \code{grouping}: Columns in \code{colData} used for the different
#'     sample groupings, accessible through \code{groupsUMI4C}.
#'     \item \code{normalized}: Logical indicating whether samples/groups are
#'      normalized or not.
#'     \item \code{region}: \linkS4class{GRanges} with the coordinates of
#'     the genomic window used for analyzing UMI4C data.
#'     \item \code{ref_umi4c}: Name of the sample or group used as reference for
#'     normalization.
#' }
#' @slot assays  Matrix where each row represents a restriction fragment site
#' and columns represent each sample or group defined in \code{grouping}.
#' After running the \code{\link{makeUMI4C}} function, it will contain the
#' following data:
#' \enumerate{
#'     \item \code{umis}: Raw number of UMIs detected by \code{\link{contactsUMI4C}}.
#'     \item \code{norm_mat}: Normalization factors for each sample/group and fragment end.
#'     \item \code{trend}: Adaptative smoothing trend of UMIs.
#'     \item \code{geo_coords}: Geometric coordinates obtained when performing
#'     the adaptative smoothing.
#'     \item \code{scale}: Scale selected for the adaptative smoothing.
#'     \item \code{sd}: Stantard deviation for the adaptative smoothing trend.
#' }
#' @slot dgram List containing the domainograms for each sample. A domainogram
#' is matrix where columns are different scales selected for merging UMI counts
#' and rows are the restriction fragments.
#' @slot groupsUMI4C List of \code{UMI4C} objects with the specific groupings.
#' @slot results List containing the results for the differential analysis ran
#' using \code{\link{fisherUMI4C}}.
#' @rdname UMI4C
#' @import methods
#' @export
.UMI4C <- setClass("UMI4C",
    slots = representation(
        dgram = "SimpleList",
        groupsUMI4C = "SimpleList",
        results = "SimpleList"
    ),
    contains = "RangedSummarizedExperiment"
)

setValidity("UMI4C", function(object) {
    TRUE
})

#' @export
#' @import SummarizedExperiment
UMI4C <- function(dgram = S4Vectors::SimpleList(),
    results = S4Vectors::SimpleList(),
    groupsUMI4C = S4Vectors::SimpleList(),
    ...) {
    se <- SummarizedExperiment(...)
    .UMI4C(se,
        dgram = dgram,
        groupsUMI4C = groupsUMI4C, 
        results = results
    )
}

#' Make UMI4C object
#'
#' The \linkS4class{UMI4C} constructor is the function \code{\link{makeUMI4C}}. By using
#'  the arguments listed below, performs the necessary steps to analyze UMI-4C
#'  data and summarize it in an object of class \linkS4class{UMI4C}.
#' @rdname UMI4C
#' @param colData Data.frame containing the information for constructing the
#' UMI4C experiment object. Needs to contain the following columns:
#' \itemize{
#'     \item sampleID. Unique identifier for the sample.
#'     \item condition. Condition for performing differential analysis. Can be
#'     control and treatment, two different cell types, etc.
#'     \item replicate. Number for identifying replicates.
#'     \item file. File as outputed by \code{umi4CatsContacts} function.
#' }
#' @param viewpoint_name Character indicating the name for the used viewpoint.
#' @param grouping Name of the column in colData used to merge the samples or
#' replicates. Set to NULL for skipping grouping. Default: "condition".
#' @param normalized Logical indicating whether UMI-4C profiles should be
#' normalized to the \code{ref_umi4c} sample/group. Default: TRUE
#' @param ref_umi4c Name of the sample or group to use as reference for
#' normalization. By default is NULL, which means it will use the sample with
#' less UMIs in the analyzed region. It should be a named vector, where the name
#' corresponds to the grouping column from \code{colData} and the value represents
#' the level to use as reference.
#' @param bait_exclusion Region around the bait (in bp) to be excluded from the
#' analysis. Default: 3000bp.
#' @param bait_expansion Number of bp upstream and downstream of the bait to use
#'  for the analysis (region centered in bait). Default: 1Mb.
#' @param scales Numeric vector containing the scales for calculating the
#' domainogram.
#' @param min_win_factor Proportion of UMIs that need to be found in a specific
#' window for adaptative trend calculation
#' @param sd Stantard deviation for adaptative trend.
#' @return It returns an object of the class \linkS4class{UMI4C}.
#' @import GenomicRanges
#' @importFrom stats as.formula
#' @seealso UMI4C-methods
#' @examples
#' # Load sample processed file paths
#' files <- list.files(system.file("extdata", "CIITA", "count",
#'     package = "UMI4Cats"
#' ),
#' pattern = "*_counts.tsv",
#' full.names = TRUE
#' )
#'
#' # Create colData including all relevant information
#' colData <- data.frame(
#'     sampleID = gsub("_counts.tsv.gz", "", basename(files)),
#'     file = files,
#'     stringsAsFactors = FALSE
#' )
#'
#' library(tidyr)
#' colData <- colData %>%
#'     separate(sampleID,
#'         into = c("condition", "replicate", "viewpoint"),
#'         remove = FALSE
#'     )
#'
#' # Load UMI-4C data and generate UMI4C object
#' umi <- makeUMI4C(
#'     colData = colData,
#'     viewpoint_name = "CIITA",
#'     grouping = "condition"
#' )
#' @export
makeUMI4C <- function(colData,
    viewpoint_name = "Unknown",
    grouping = "condition",
    normalized = TRUE,
    ref_umi4c = NULL,
    bait_exclusion = 3e3,
    bait_expansion = 1e6,
    scales = 5:150,
    min_win_factor = 0.02,
    sd = 2) {
    if (!("condition" %in% names(colData))) {
          stop("colData must contain 'condition'")
      }
    if (!("replicate" %in% names(colData))) {
          stop("colData must contain 'replicate'")
    }
    if (!("sampleID" %in% names(colData))) {
      stop("colData must contain 'sampleID'")
    }
    if (!("file" %in% names(colData))) {
          stop("colData must contain 'file'")
    }
    if (length(grouping)>1) {
      stop("Use only one varible for grouping. You can latter add more groupings using the addGrouping() function.")
    }
    if (!is.null(grouping)) {
      if(!(grouping %in% colnames(colData))) {
        stop("Grouping variable not found among colnames(colData).")
      }
    }

    colData$sampleID <- gsub(".", "_", colData$sampleID, fixed = TRUE)

    ## Load UMI4C matrices
    mats <- lapply(as.character(colData$file),
        utils::read.delim,
        header = TRUE,
        stringsAsFactors = FALSE
    )
    names(mats) <- colData$sampleID

    nrows <- vapply(mats, nrow, FUN.VALUE=integer(1))
    if (length(unique(nrows)) != 1) stop("Number of rows for the supplied files is different. Please check again your methods.")

    pos <- lapply(mats, function(x) paste0(x[, "chr_contact"], ":", x[, "start_contact"], "-", x[, "end_contact"]))
    if (length(unique(pos)) != 1) stop("Fragment end coordinates for your files are different. Please check again your methods.")

    max <- lapply(mats, function(x) max(x[, "UMIs"], na.rm = TRUE))
    is_zero <- names(max)[max == 0]
    if (length(is_zero) > 0) {
        message(
            "Warning:\n",
            "Your samples ", paste0(is_zero, collapse = " "),
            " don't have any UMIs for the given fragment ends. ",
            "Check again your analysis and experiment quality."
        )
    }

    ## Obtain bait information
    baits <- lapply(
        mats,
        function(x) {
            unique(GRanges(
                seqnames = unique(x[, "chr_bait"]),
                ranges = IRanges(
                    start = x[, "start_bait"],
                    end = x[, "end_bait"]
                )
            ))
        }
    )
    baits <- unlist(GRangesList(baits))
    bait <- unique(baits)
    if (length(bait) > 1) {
          stop("Bait position for the supplied samples differ. Please check again your methods.")
      }

    bait$name <- viewpoint_name
    rownames(bait) <- NULL

    ## Create assay matrix with raw UMIs
    umis <- do.call(rbind, mats)
    umis <- umis[, !grepl("bait", colnames(umis))]
    colnames(umis) <- c("chr_contact", "start_contact", "end_contact", "UMIs")
    umis$sampleID <- unlist(lapply(
        strsplit(rownames(umis), ".", fixed = TRUE),
        function(x) x[1]
    ))
    umis.d <- reshape2::dcast(umis,
        chr_contact + start_contact + end_contact ~ sampleID,
        value.var = "UMIs"
    )

    umis.d <- umis.d[order(umis.d$chr_contact, umis.d$start_contact), ]
    umis.d$id_contact <- paste0("frag_", seq_len(nrow(umis.d)))

    # Create row_ranges
    rowRanges <- GRanges(
        seqnames = umis.d$chr_contact,
        ranges = IRanges(
            start = umis.d$start_contact,
            end = umis.d$end_contact
        )
    )
    rowRanges$id_contact <- umis.d$id_contact

    # Create assay matrix
    assay <- as.matrix(umis.d[,-c(seq_len(3), ncol(umis.d))])
    rownames(assay) <- umis.d$id_contact

    ## Create summarizedExperiment
    umi4c <- UMI4C(
        colData = colData,
        rowRanges = rowRanges,
        metadata = list(
            bait = bait,
            scales = scales,
            min_win_factor = min_win_factor,
            normalized = normalized,
            grouping = NULL
        ),
        assays = SimpleList(umi = assay)
    )

    ## Remove region around bait
    bait_exp <- GenomicRanges::resize(metadata(umi4c)$bait,
        fix = "center",
        width = bait_exclusion
    )
    umi4c <- subsetByOverlaps(umi4c, bait_exp, invert = TRUE)

    if (any(colSums(assay(umi4c)) == 0)) {
        stop("The number of UMIs at least for one sample are 0. Try reducing your
         bait_exclusion value.")
    }

    ## Remove regions outside scope
    region <- GenomicRanges::resize(metadata(umi4c)$bait,
        fix = "center",
        width = bait_expansion
    )
    umi4c <- subsetByOverlaps(umi4c, region)
    metadata(umi4c)$region <- region

    ## Divide upstream and downstream coordinates
    rowRanges(umi4c)$position <- NA
    rowRanges(umi4c)$position[start(rowRanges(umi4c)) < start(metadata(umi4c)$bait)] <- "upstream"
    rowRanges(umi4c)$position[start(rowRanges(umi4c)) >= start(metadata(umi4c)$bait)] <- "downstream"
    
    ##### PROCESSING FOR PLOTTING
    ## Add stats using sampleIDs -------------
    ref <- metadata(umi4c)$ref_umi4c
    
    if (is.null(ref) | !("sampleID" %in% names(ref))) {
      # Get sample with less UMIs if no ref present
      metadata(umi4c)$ref_umi4c <- colnames(assay(umi4c))[which(colSums(assay(umi4c)) == min(colSums(assay(umi4c))))]
    } else {
      # Use value from named list
      metadata(umi4c)$ref_umi4c <- refs[grouping]
    }
    
    # Get normalization matrix
    umi4c <- getNormalizationMatrix(umi4c)
    
    ## Calculate domainograms
    umi4c <- calculateDomainogram(umi4c, scales = scales, normalized = normalized)
    
    ## Calculate adaptative trend
    umi4c <- calculateAdaptativeTrend(umi4c, sd = sd, normalized = normalized)

    ## Use groupings -----
    if (!is.null(grouping)) {
      umi4c <- addGrouping(umi4c, grouping=grouping, normalized=normalized, scales=scales, sd=sd)
    }
    
    return(umi4c)
}