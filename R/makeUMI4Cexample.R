#' Make UMI4C example object
#'
#' @param ... Other arguments to be passed to the \code{makeUMI4C} function.
#' @return Creates a UMI4C object to use as an example.
#' @examples
#' umi4c <- makeUMI4Cexample()
#' @importFrom rlang .data
#' @export
makeUMI4Cexample <- function(...) {
  # Load sample processed file paths
  files <- list.files(system.file("extdata", "SOCS1", "count",
                                  package="UMI4Cats"),
                      pattern="*_counts.tsv.gz",
                      full.names=TRUE)

  # Create colData including all relevant information
  colData <- data.frame(sampleID = gsub("_counts.tsv.gz", "", basename(files)),
                        file = files,
                        stringsAsFactors=FALSE)

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

#' Download UMI4Cats example datasets
#'
#' Downloads the required UMI4Cats example datasets.
#' @param output_dir Output directory for the datasets.
#' @param file_dir Path to the compressed IRB database file.
#' @return It creates the \code{output_dir} with the example UMI-4C files used by the vignette.
#' @examples
#' path <- downloadUMI4CexampleData()
#' @importFrom utils download.file untar
#' @export
downloadUMI4CexampleData <- function(output_dir="./",
                                     file_dir="http://gattaca.imppc.org/genome_browser/lplab/UMI4Cats_data.tar.gz") {
  tf <- tempfile()
  message("Will begin downloading datasets to ", tf)
  ret <- download.file(file_dir, tf, mode='wb')

  if (ret != 0) stop("Couldn't download file from ", file_dir)

  untar(tf, exdir=output_dir, verbose=TRUE)
  message("Done writing UMI4C example files to ", output_dir)

  return(invisible(paste0(output_dir, "UMI4Cats_data/")))
}
