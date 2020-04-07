#' UMI4Cats: A package for analyzing UMI-4C chromatin contact data
#'
#' The UMI4Cats package provides functions for the pre-processing, analysis and
#' visualization of UMI-4C chromatin contact data.
#'
#' @section Pre-processing:
#' The pre-processing functions ...
#'
#' @section Analysis:
#' The analysis functions ...
#'
#' @section Visualization:
#' The plotting functions ...
#'
#' @docType package
#' @name UMI4Cats
NULL

utils::globalVariables(c("factors", "scales", "value",
                         "hg19_gene_annoation_ensemblv75", "stepping",
                         "gene_name", "geo_coord", "grouping_var",
                         "relative_position", "sample_id", "al_mapped",
                         "al_unmapped", "variable",
                         "queryHits", "subjectHits", "subsetByOverlaps",
                         "assay", "metadata", "rowRanges", "assays", "assays<-",
                         "assay", "colData", "IRanges", "SimpleList",
                         "metadata<-", "rowRanges<-", "UMIs"))
