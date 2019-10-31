#' @export
setMethod("dgram", "UMI4C",
          function(x, ...) {
            out <- x@dgram

            out
          })

#' @exportMethod "dgram<-"
setReplaceMethod("dgram",
                 signature(object="UMI4C", value="SimpleList"),
                 function(object, value) {

                   object@dgram <- value
                   return(object)
                 })

#' @export
setMethod("bait", "UMI4C",
          function(x, ...) {
            out <- x@metadata$bait

            out
          })

#' @export
setMethod("trend", "UMI4C",
          function(x, ...) {
            ## Construct trend df using geo_coords and trend
            trend_df <- data.frame(geo_coord = as.vector(assays(x)$geo_coord),
                                   trend = as.vector(assays(x)$trend),
                                   sd = as.vector(assays(x)$sd),
                                   scale = as.vector(assays(x)$scale),
                                   id_contact = rep(rowRanges(x)$id_contact,
                                                    ncol(x)))
            trend_df$sample <- rep(colnames(assay(x)),
                                   each=nrow(x))

            trend_df <- trend_df[!is.na(trend_df$trend),]

            trend_df <- dplyr::left_join(trend_df,
                                         data.frame(colData(x)),
                                         by=c(sample="sampleID"))
            trend_df
          })

#' @export
setMethod("results", "UMI4C",
          function(x, format="GRanges", counts=FALSE, ...) {
            res_list <- x@results

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

            results
          })
