#' Demultiplex FASTQ files using fastq-multx
#'
#' Demultiplex FASTQ files containng different bait information
#' @param barcodes Dataframe with "name of sample" and "barcode" for every
#' sample to demultiplex
#' @param fastq Fastq to demultiplex containing mate 1s. Different pairs should
#' be named "_R1" or "_R2". Allowed formats: _R1.fastq.gz, _R1.fq.gz, _R1.fastq
#' or _R1.fq.
#' @param numb_reads Number of lines from the FastQ file to load in each loop.
#' If having memory size problems, change it to a smaller number. Default=10e10.
#' @param out_path Path where to save the demultiplex output. Defaults to a path
#'  named \code{raw_fastq} in your working
#' directory.
#' @return Paired-end FastQ files compressed.
#' @examples
#' \dontrun{
#' path <- downloadUMI4CexampleData(use_sample=TRUE)
#' fastq <- file.path(path, "SOCS1", "fastq", "sub_ctrl_hi19_SOCS1_R1.fastq.gz")
#' barcodes <- data.frame(sample=c("SOCS1"),
#'                        barcode=c("CCCAAATCGCCCAGACCAG"))
#'
#' demultiplexFastq(barcodes=barcodes,
#'                  fastq=fastq,
#'                  out_path = path)
#'}
#'
#'@export
demultiplexFastq <- function(barcodes,
                             fastq,
                             out_path="raw_fastq",
                             numb_reads=10e10){
  fq_R1 <- fastq
  fq_R2 <- gsub('_R1.', '_R2.', fq_R1)

  extension_fastq <- utils::tail(unlist(strsplit(fq_R1, "_")), n=1)

  if(!extension_fastq %in% c('R1.fastq.gz', 'R1.fq.gz', 'R1.fastq', 'R1.fq')){
    stop(paste("FASTQ file should have one of the following extensions:
               _R1.fastq.gz, _R1.fq.gz, _R1.fastq or _R1.fq"))
  }

  if(length(fq_R2) == 0) stop(paste("Files should be paired-end type"))

  if(!is(barcodes, "data.frame")) stop(paste("Barcodes should be a Dataframewith name of sample and barcode for every sample to demultiplex"))

  message(paste("Starting demultiplex using:\n",
                "> Barcodes:\n", paste(utils::capture.output(print(barcodes)),
                                       collapse = "\n"), "\n\n",
                "> File R1:", fq_R1, "\n",
                "> File R2:", fq_R2, "\n",
                "> Output path:", out_path))



  # Initialize global variables
  total_reads <- 0
  specific_reads <- 0


  for (i in seq_len(nrow(barcodes))){

    stream1 <- ShortRead::FastqStreamer(fq_R1)
    stream2 <- ShortRead::FastqStreamer(fq_R2)

    # generate barcode
    barcode <- barcodes$barcode[i]

    repeat {
      reads_fqR1 <- ShortRead::yield(stream1, n = numb_reads)
      reads_fqR2 <- ShortRead::yield(stream2, n = numb_reads)

      if (length(reads_fqR1)==0) break
      if (length(reads_fqR1)!=length(reads_fqR2)) stop("Different number of reads in R1 vs R2")

      total_reads <- total_reads + length(reads_fqR1) # Save total reads

      # for cases when the bait is to far from restriction enzyme
      if (nchar(as.character(barcode)) > unique(ShortRead::width(reads_fqR1))){
        barcode = substr(barcode, 1, unique(width(reads_fqR1)))
      }

      # filter reads that not present barcode
      barcode_reads_fqR1 <- reads_fqR1[grepl(barcode, ShortRead::sread(reads_fqR1))]
      barcode_reads_fqR2 <- reads_fqR2[grepl(barcode, ShortRead::sread(reads_fqR1))]

      specific_reads <- specific_reads + length(barcode_reads_fqR1) # Save specific reads

      # write output fastq files

      output_fastq <- file.path(out_path, barcodes$sample[i])

      ShortRead::writeFastq(barcode_reads_fqR1,
                            paste0(output_fastq, '_R1.fq.gz'),
                            mode="a")

      ShortRead::writeFastq(barcode_reads_fqR2,
                            paste0(output_fastq, '_R2.fq.gz'),
                            mode="a")
    }

    # Construct stats data.frame
    stats <- data.frame(sample_id=barcodes$sample[i],
                        total_reads=total_reads,
                        specific_reads=specific_reads,
                        stringsAsFactors = FALSE)

    # create stats file and save
    stats <- do.call(rbind, stats)
    utils::write.table(stats,
                       file = file.path(out_path,
                                        paste0(barcodes$sample[i],
                                               "_umi4cats_demultiplexFastq_stats.txt")),
                       row.names = FALSE,
                       sep="\t",
                       quote=FALSE)

    message("Finished demultiplex sample ", barcodes$sample[i])
  }

  on.exit(close(stream1))
  on.exit(close(stream2), add =TRUE)
  }


