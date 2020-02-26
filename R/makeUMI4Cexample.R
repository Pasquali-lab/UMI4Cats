#' Make UMI4C example object
#'
#' @param ... Other arguments to be passed to the \code{makeUMI4C} function.
#' @examples
#' umi4c <- makeUMI4Cexample()
#' @importFrom rlang .data
#' @return Creates a UMI4C object to use as an example.
#' @export
makeUMI4Cexample <- function(...) {
  # Load sample processed file paths
  files <- list.files(system.file("extdata", "SOCS1", "count",
                                  package="UMI4Cats"),
                      pattern="*_counts.tsv",
                      full.names=TRUE)

  # Create colData including all relevant information
  colData <- data.frame(sampleID = gsub("_counts.tsv", "", basename(files)),
                        file = files,
                        stringsAsFactors=F)

  colData <- colData %>%
    tidyr::separate(.data$sampleID,
             into=c("condition", "replicate", "viewpoint"),
             remove=FALSE)

  # Load UMI-4C data and generate UMI4C object
  umi <- makeUMI4C(colData=colData,
                   viewpoint_name="SOCS1",
                   ...)

  return(umi)
}
