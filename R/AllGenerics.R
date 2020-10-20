#' @rdname UMI4C-methods
#' @export
setGeneric("dgram", function(object) standardGeneric("dgram"))

#' @rdname UMI4C-methods
#' @export
setGeneric("dgram<-", function(object, value) standardGeneric("dgram<-"))

#' @rdname UMI4C-methods
#' @export
setGeneric("groupsUMI4C", function(object, value) standardGeneric("groupsUMI4C"))

#' @rdname UMI4C-methods
#' @export
setGeneric("groupsUMI4C<-", function(object, value) standardGeneric("groupsUMI4C<-"))

#' @rdname UMI4C-methods
#' @export
setGeneric("bait", function(object) standardGeneric("bait"))

#' @rdname UMI4C-methods
#' @export
setGeneric("trend", function(object) standardGeneric("trend"))

#' @rdname UMI4C-methods
#' @export
setGeneric(
    "resultsUMI4C",
    function(object,
    format = "GRanges",
    counts = TRUE,
    ordered = FALSE) {
        standardGeneric("resultsUMI4C")
    }
)
