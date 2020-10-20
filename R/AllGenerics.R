#' @rdname UMI4C-methods
#' @export
setGeneric("dgram", function(object) standardGeneric("dgram"))

#' @rdname UMI4C-methods
#' @export
setGeneric("dgram<-", function(object, value) standardGeneric("dgram<-"))

#' @rdname UMI4C-methods
#' @export
setGeneric("bait", function(object) standardGeneric("bait"))

#' @rdname UMI4C-methods
#' @export
setGeneric("trend", function(object) standardGeneric("trend"))

#' @rdname UMI4C-methods
#' @export
setGeneric(
    "results",
    function(object,
    format = "GRanges",
    counts = TRUE,
    ordered = FALSE) {
        standardGeneric("results")
    }
)
