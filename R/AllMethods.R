#' @title UMI4C class methods
#' @docType methods
#'
#' @description This page contains a summary of the different methods used to
#' access the information contained inside the \code{UMI4C-class} object. See
#' the details section for more information on the different accessors.
#' @name UMI4C-methods
#'
#' @param object a \code{UMI4C-class} object.
#'
#' @return There are several accessors to easily retrive information from a
#' \code{UMI4C-class} object:
#' \itemize{
#'   \item \code{dgram}: Returns a named list with the output domainograms for
#'   each sample.
#'   \item \code{bait}: Returns a \linkS4class{GRanges} object with the position
#'    of the bait.
#'   \item \code{trend}: Returns a data.frame in long format with the values of
#'   the adapative smoothen trend.
#'   \item \code{resultsUMI4C}: Returns a \linkS4class{GRanges} or data.frame with
#'   the results of the differential analysis.
#' }
#' @examples
#' # Access the different information inside the UMI4C object
#' data("ex_ciita_umi4c")
#' ex_ciita_umi4c <- addGrouping(ex_ciita_umi4c, grouping="condition")
#' 
#' dgram(ex_ciita_umi4c)
#' bait(ex_ciita_umi4c)
#' head(trend(ex_ciita_umi4c))
#'
#' # Perform differential test
#' enh <- GRanges(c("chr16:10925006-10928900", "chr16:11102721-11103700"))
#' umi_dif <- fisherUMI4C(ex_ciita_umi4c, query_regions = enh, 
#'                        filter_low = 20, resize = 5e3)
#' resultsUMI4C(umi_dif)
#' @seealso UMI4C, UMI4C-class
NULL

#' @name dgram
#' @param object a \code{UMI4C-class} object.
#' @rdname UMI4C-methods
#' @aliases dgram,UMI4C-method
#'
#' @export
setMethod(
    "dgram", "UMI4C",
    function(object) {
        out <- object@dgram

        out
    }
)

#' @name dgram
#' @param object a \code{UMI4C-class} object.
#' @param value Alternative list of dgrams to replace the current slot.
#' @rdname UMI4C-methods
#' @aliases dgram<-,UMI4C-method
#' @exportMethod "dgram<-"
setReplaceMethod(
    "dgram",
    "UMI4C",
    function(object, value) {
        object@dgram <- value
        return(object)
    }
)

#' @name groupsUMI4C
#' @param object a \code{UMI4C-class} object.
#' @rdname UMI4C-methods
#' @aliases groupsUMI4C,UMI4C-method
#'
#' @export
setMethod(
    "groupsUMI4C", "UMI4C",
    function(object) {
        out <- object@groupsUMI4C
        
        out
    }
)

#' @name groupsUMI4C
#' @param object a \code{UMI4C-class} object.
#' @param value Alternative list of dgrams to replace the current slot.
#' @rdname UMI4C-methods
#' @aliases groupsUMI4C<-,UMI4C-method
#' @exportMethod "groupsUMI4C<-"
setReplaceMethod(
    "groupsUMI4C",
    "UMI4C",
    function(object, value) {
        object@groupsUMI4C <- value
        return(object)
    }
)

#' @name bait
#' @rdname UMI4C-methods
#' @aliases bait,UMI4C-method
#' @export
setMethod(
    "bait", "UMI4C",
    function(object) {
        out <- object@metadata$bait

        out
    }
)

#' @name trend
#' @param object a \code{UMI4C-class} object.
#' @rdname UMI4C-methods
#' @aliases trend,UMI4C-method
#' @export
setMethod(
    "trend", "UMI4C",
    function(object) {
        ## Construct trend df using geo_coords and trend
        trend_df <- data.frame(
            geo_coord = as.vector(assays(object)$geo_coord),
            trend = as.vector(assays(object)$trend),
            sd = as.vector(assays(object)$sd),
            scale = as.vector(assays(object)$scale),
            id_contact = rep(
                rowRanges(object)$id_contact,
                ncol(object)
            )
        )
        trend_df$sample <- rep(colnames(assay(object)),
            each = nrow(object)
        )

        trend_df <- trend_df[!is.na(trend_df$trend), ]
        
        group <- metadata(object)$grouping
        if (is.null(group)) group <-  "sampleID"

        trend_df <- dplyr::left_join(trend_df,
            data.frame(colData(object)),
            by = c(sample = group)
        )
        trend_df
    }
)

#' @name resultsUMI4C
#' @param object a \code{UMI4C-class} object.
#' @param format Either "GRanges" (default) or "data.frame", indicating the
#' format output of the results.
#' @param counts Logical indicating whether counts for the different region
#' should be provided. Default: FALSE.
#' @param ordered Logical indicating whether to sort output by significance
#' (adjusted p-value). Default: FALSE.
#' @rdname UMI4C-methods
#' @aliases resultsUMI4C,UMI4C-method
#' @export
setMethod(
    "resultsUMI4C", "UMI4C",
    function(object, format = "GRanges", counts = FALSE, ordered = FALSE) {
        res_list <- object@results

        results <- res_list$query

        if (counts) {
            mcols(results) <- dplyr::left_join(data.frame(mcols(results)),
                res_list$counts,
                by = c(id = "query_id")
            )
        }

        mcols(results) <- dplyr::left_join(data.frame(mcols(results)),
            res_list$results,
            by = c(id = "query_id")
        )

        names(results) <- NULL

        if (format == "data.frame") results <- data.frame(results)[, -c(4:5)]

        if (ordered) {
            results <- results[order(results$padj,
                decreasing = FALSE
            ), ]
        }

        results
    }
)
