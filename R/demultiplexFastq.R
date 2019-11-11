#' Demultiplex FASTQ files using fastq-multx
#'
#' Demultiplex FASTQ files containng different bait information using \link[=https://github.com/brwnj/fastq-multx]{fastq-multx} software
#' @param barcodes Dataframe with "name of sample" and "barcode" for every sample to demultiplex
#' @param fastq Name of the fastq to demultiplex. When paired-end, you only need to provide the name for the "_R1" file. Different pairs should be named "_R1" or "_R2".
#' @param type Character indicating wether fastq files are in "single" or "paired" end format. Default is \code{"paired"} format.
#' @param out_path Path where to save the demultiplex output. Defaults to a path named \code{raw_fastq} in your working
#' directory.
#' @examples
#' \dontrun{
#' barcodes <- data.frame(sample=c("KLK3", "KLK6", "KLK7", "KLK10"),
#'                        barcode=c("ATGGTCTGGGCGCTGTCTTG",
#'                                  "TATTCTTCCTCAGCCCACATCTT",
#'                                  "GGATGAAGATTTTGGAGCCCAGC",
#'                                  "GGGCGGGGATTGAACGC"))
#'
#' demultiplexFastq(barcodes=barcodes,
#'                  fastq="~/samples/sample_1_R1.fastq.gz",
#'                  type="paired",
#'}
#'
#'@export

demultiplexFastq <- function(barcodes,
                             fastq,
                             type="paired",
                             out_path="raw_fastq"){

  barcode_file <- file.path(out_path, "barcodes.txt")

  write.table(barcodes,
              barcode_file,
              col.names = F, row.names = F, quote = F,
              sep = "\t")

  # TODO: This won't work because binary is not include in the package
  # fastqMultx <- system.file("bin/fastq-multx/fastq-multx",
  #                           package = "UMI4Cats")

  fastq_multx <- "fastq_multx" # workaround for now

  name_file <- gsub("\\..*$", "", basename(fastq))

  message(paste("Starting demultiplex using:\n",
          "> Barcodes:\n", paste(capture.output(print(barcodes)), collapse = "\n"), "\n\n",
          "> File:", fastq, "\n",
          "> Type:", type, "\n",
          "> Output path:", out_path))


  if (type == "paired") { # file paired end
    file_R1 <- fastq[grep("_R1", fastq)]
    file_R2 <- unique(gsub("_R1", "_R2", fastq))

    out_R1 <- file.path(out_path, "%_R1.fq.gz")
    out_R2 <- file.path(out_path, "%_R2.fq.gz")

    system(paste(fastqMultx,
                 "-x -m 0",
                 "-b", barcode_file,
                 file_R1, file_R2,
                 "-o", out_R1, "-o", out_R2))

    } else if (type == "single") { # file paired end
      out_file <- file.path(out_path, "%.fq.gz")

      system(paste(fastqMultx,
                   "-x -m 0",
                   "-b", barcode_file,
                   fastq,
                   out_file))
    }

  unlink(barcode_file)

  message("Finished demultiplex")
}
