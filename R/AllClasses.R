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
#'     \item \code{grouping}: Column in \code{colData} used to group the samples.
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
#' @slot results List containing the results for the differential analysis ran
#' using \code{\link{fisherUMI4C}}.
#' @rdname UMI4C
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @export
.UMI4C <- setClass("UMI4C",
    slots = representation(
        dgram = "SimpleList",
        results = "SimpleList"
    ),
    contains = "RangedSummarizedExperiment"
)

setValidity("UMI4C", function(object) {
    TRUE
})

#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
UMI4C <- function(dgram = S4Vectors::SimpleList(),
    results = S4Vectors::SimpleList(),
    ...) {
    se <- SummarizedExperiment(...)
    .UMI4C(se,
        dgram = dgram,
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
#' replicates. Default: "condition".
#' @param normalized Logical indicating whether UMI-4C profiles should be
#' normalized to the \code{ref_umi4c} sample/group. Default: TRUE
#' @param ref_umi4c Name of the sample or group to use as reference for
#' normalization. By default is NULL, which means it will use the sample with
#' less UMIs in the analyzed region. It should be one of the values from the
#' column used as \code{grouping}.
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
#' @importFrom SummarizedExperiment SummarizedExperiment
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
    if (!("file" %in% names(colData))) {
          stop("colData must contain 'file'")
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

    pos <- lapply(mats, function(x) paste0(x[, 3], ":", x[, 4]))
    if (length(unique(pos)) != 1) stop("Fragment end coordinates for your files are different. Please check again your methods.")

    max <- lapply(mats, function(x) max(x[, 5], na.rm = TRUE))
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
                seqnames = unique(x[, 1]),
                ranges = IRanges(
                    start = x[, 2],
                    end = x[, 2]
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
    umis <- umis[, -c(1, 2)]
    colnames(umis) <- c("chr_contact", "pos_contact", "UMIs")
    umis$sampleID <- unlist(lapply(
        strsplit(rownames(umis), ".", fixed = TRUE),
        function(x) x[1]
    ))
    umis.d <- reshape2::dcast(umis,
        chr_contact + pos_contact ~ sampleID,
        value.var = "UMIs"
    )

    umis.d <- umis.d[order(umis.d$chr_contact, umis.d$pos_contact), ]
    umis.d$id_contact <- paste0("frag_", seq_len(nrow(umis.d)))

    # Create row_ranges
    rowRanges <- GRanges(
        seqnames = umis.d$chr_contact,
        ranges = IRanges(
            start = umis.d$pos_contact,
            end = umis.d$pos_contact
        )
    )
    rowRanges$id_contact <- umis.d$id_contact

    # Create assay matrix
    if (grouping %in% colnames(colData)) { ## Create assays for grouping variables

        ## Sum UMI4C from replicates
        assay_all <- as.matrix(umis.d[, -c(1, 2, ncol(umis.d)), drop = FALSE], labels = TRUE)
        rownames(assay_all) <- rowRanges$id_contact

        assay_m <- reshape2::melt(assay_all)
        colnames(assay_m) <- c("rowname", "sampleID", "UMIs")
        assay_m <- suppressWarnings(dplyr::left_join(
            assay_m,
            colData[, unique(c("sampleID", grouping)), drop = FALSE],
            by = "sampleID"
        ))
        assay_df <- assay_m %>%
            dplyr::group_by_at(c("rowname", grouping)) %>%
            dplyr::summarise(UMIs = sum(UMIs, na.rm = TRUE)) %>%
            reshape2::dcast(stats::as.formula(paste0("rowname~", grouping)), value.var = "UMIs")

        assay <- as.matrix(assay_df[, -which(colnames(assay_df) == "rowname")], )
        rownames(assay) <- assay_df$rowname
        colnames(assay) <- colnames(assay_df)[-1]

        ## Summarize colData
        colData <- colData %>%
            dplyr::group_by_at(grouping) %>%
            dplyr::summarise_all(paste0, collapse = ", ")
    }

    ## Create summarizedExperiment
    umi4c <- UMI4C(
        colData = colData,
        rowRanges = rowRanges,
        metadata = list(
            bait = bait,
            scales = scales,
            min_win_factor = min_win_factor,
            grouping = grouping,
            normalized = normalized
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

    ## Get normalization matrix
    if (is.null(ref_umi4c)) {
        metadata(umi4c)$ref_umi4c <- colnames(assay(umi4c))[which(colSums(assay(umi4c)) == min(colSums(assay(umi4c))))]
    } else {
        ref_umi4c <- gsub(".", "_", ref_umi4c, fixed = TRUE)

        if (!(ref_umi4c %in% dplyr::pull(colData, grouping))) {
            stop("The name of the sample in ref_umi4c does not correspond to a sample name in colData.")
        }

        metadata(umi4c)$ref_umi4c <- ref_umi4c
    }


    umi4c <- getNormalizationMatrix(umi4c)

    ## Calculate domainograms
    umi4c <- calculateDomainogram(umi4c,
        scales = scales,
        normalized = normalized
    )

    ## Calculate adaptative trend
    umi4c <- calculateAdaptativeTrend(umi4c,
        sd = sd,
        normalized = normalized
    )

    return(umi4c)
}
