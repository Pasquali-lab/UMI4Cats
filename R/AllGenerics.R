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
setGeneric("results", function(object, format="GRanges", counts=TRUE) standardGeneric("results"))


#' @importMethodsFrom SummarizedExperiment assay
#' @export
setGeneric("assay")

#' @importMethodsFrom SummarizedExperiment assays
#' @export
setGeneric("assays")

#' @importMethodsFrom SummarizedExperiment assays<-
#' @export
setGeneric("assays<-")

#' @importMethodsFrom SummarizedExperiment colData
#' @export
setGeneric("colData")

#' @importMethodsFrom SummarizedExperiment colData<-
#' @export
setGeneric("colData<-")

#' @importMethodsFrom S4Vectors metadata
#' @export
setGeneric("metadata")

#' @importMethodsFrom S4Vectors metadata<-
#' @export
setGeneric("metadata<-")

#' @importMethodsFrom SummarizedExperiment rowRanges
#' @export
setGeneric("rowRanges")

#' @importMethodsFrom SummarizedExperiment rowRanges<-
#' @export
setGeneric("rowRanges<-")

#' @importFrom S4Vectors queryHits
#' @export
setGeneric("queryHits")

#' @importFrom S4Vectors subjectHits
#' @export
setGeneric("subjectHits")

#' @importFrom IRanges subsetByOverlaps
#' @export
setGeneric("subsetByOverlaps")

#' @importFrom IRanges IRanges
#' @export
setGeneric("IRanges")

#' @importFrom S4Vectors SimpleList
#' @export
setGeneric("SimpleList")

#' @importMethodsFrom  S4Vectors elementNROWS
#' @export
setGeneric("elementNROWS")

