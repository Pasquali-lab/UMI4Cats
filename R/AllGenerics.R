#' @export
setGeneric("dgram", function(x, ...) standardGeneric("dgram"))

#' @export
setGeneric("dgram<-", function(object, value, ...) standardGeneric("dgram<-"))

#' @export
setGeneric("bait", function(x, ...) standardGeneric("bait"))

#' @export
setGeneric("trend", function(x, ...) standardGeneric("trend"))

#' @export
setGeneric("results", function(x, format="GRanges", counts=TRUE, ...) standardGeneric("results"))
