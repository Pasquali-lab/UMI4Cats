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
#'   \item \code{results}: Returns a \linkS4class{GRanges} or data.frame with
#'   the results of the differential analysis.
#' }
#' @examples
#' # Create example UMI4C-class object with helper function
#' umi <- makeUMI4Cexample()
#'
#' # Access the different information inside the UMI4C object
#' dgram(umi)
#' bait(umi)
#' trend(umi)
#'
#' # Perform differential test
#' umi <- fisherUMI4C(umi)
#' results(umi)
#'
#' @seealso UMI4C, UMI4C-class
NULL

#' @name dgram
#' @param object a \code{UMI4C-class} object.
#' @rdname UMI4C-methods
#' @aliases dgram,UMI4C-method
#'
#' @export
setMethod("dgram", "UMI4C",
          function(object) {
            out <- object@dgram

            out
          })

#' @name dgram
#' @param object a \code{UMI4C-class} object.
#' @param value Alternative list of dgrams to replace the current slot.
#' @rdname UMI4C-methods
#' @aliases dgram<-,UMI4C-method
#' @exportMethod "dgram<-"
setReplaceMethod("dgram",
                 "UMI4C",
                 function(object, value) {

                   object@dgram <- value
                   return(object)
                 })

#' @name bait
#' @rdname UMI4C-methods
#' @aliases bait,UMI4C-method
#' @export
setMethod("bait", "UMI4C",
          function(object) {
            out <- object@metadata$bait

            out
          })

#' @name trend
#' @param object a \code{UMI4C-class} object.
#' @rdname UMI4C-methods
#' @aliases trend,UMI4C-method
#' @export
setMethod("trend", "UMI4C",
          function(object) {
            ## Construct trend df using geo_coords and trend
            trend_df <- data.frame(geo_coord = as.vector(assays(object)$geo_coord),
                                   trend = as.vector(assays(object)$trend),
                                   sd = as.vector(assays(object)$sd),
                                   scale = as.vector(assays(object)$scale),
                                   id_contact = rep(rowRanges(object)$id_contact,
                                                    ncol(object)))
            trend_df$sample <- rep(colnames(assay(object)),
                                   each=nrow(object))

            trend_df <- trend_df[!is.na(trend_df$trend),]

            trend_df <- dplyr::left_join(trend_df,
                                         data.frame(colData(object)),
                                         by=c(sample=metadata(object)$grouping))
            trend_df
          })

#' @name results
#' @param object a \code{UMI4C-class} object.
#' @param format Either "GRanges" (default) or "data.frame", indicating the
#' format output of the results.
#' @param counts Logical indicating whether counts for the different region
#' should be provided. Default: FALSE.
#' @param ordered Logical indicating whether to sort output by significance
#' (adjusted p-value). Default: FALSE.
#' @rdname UMI4C-methods
#' @aliases results,UMI4C-method
#' @export
setMethod("results", "UMI4C",
          function(object, format="GRanges", counts=FALSE, ordered=FALSE) {
            res_list <- object@results

            results <- res_list$query

            if (counts) {
              mcols(results) <- dplyr::left_join(data.frame(mcols(results)),
                                                 res_list$counts,
                                                 by=c(id="query_id"))
            }

            mcols(results) <- dplyr::left_join(data.frame(mcols(results)),
                                               res_list$results,
                                               by=c(id="query_id"))

            names(results) <- NULL

            if (format=="data.frame") results <- data.frame(results)[,-c(4:5)]

            if (ordered) results <- results[order(results$padj,
                                                  decreasing=FALSE)]

            results
          })



