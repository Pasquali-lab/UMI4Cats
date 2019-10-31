#' Demultiplex using fastq-multx
#'
#'@description
#' Demultiplex samples using fastq-multx software
#'
#'@usage
#'fastqMultxR(dfBarcodes, fastqFile, type, outPath)
#'
#'@param dfBarcodes Dataframe with "name of sample" and "barcode" for every sample to demultiplex
#'@param fastqFile Fastq to demultiplex.
#'In paired-end file used should be R1 and following the format "file_R1.*"
#'@param type "single" or "paired" format
#'@param outPath Path where to save the demultiplex output
#'
#'@examples
#'\dontrun{
#'barcodeM <- c('KLK3',	'ATGGTCTGGGCGCTGTCTTG',
#'              'KLK6',	'TATTCTTCCTCAGCCCACATCTT',
#'              'KLK7',	'GGATGAAGATTTTGGAGCCCAGC',
#'              'KLK10',	'GGGCGGGGATTGAACGC')
#'
#'barcodeM <- t(matrix(barcodeM, nrow = 2))
#'dfBarcodes <- data.frame(barcodeM)
#'
#'colnames(dfBarcodes) <- c('sample', 'barcode')
#'fastqFile <- '/home/labs/lplab/msubirana/Documents/IGTP/kallikreins/raw/fastq/umi4c/CTTGTA_R1.fastq.gz'
#'type <- 'paired'
#'outPath <- '/home/labs/lplab/msubirana/Documents/IGTP/kallikreins/raw/fastq'
#'
#'fastqMultxR(dfBarcodes,
#'            fastqFile,
#'            type,
#'            outPath)
#'}
#'
#'@export


fastqMultxR <- function(dfBarcodes,
                   fastqFile,
                   type,
                   outPath){

  barcodeFile <- file.path(outPath, "barcodes.txt")

  write.table(dfBarcodes,
              barcodeFile,
              col.names = F, row.names = F, quote = F,
              sep = "\t")

  fastqMultx <- system.file("bin/fastq-multx/fastq-multx", package = "UMI4Cats")

  nameFile <- gsub("\\..*$", "", basename(fastqFile))

  message(paste('Starting demultiplex using:\n',
          'Barcodes:', "\n\n", paste(capture.output(print(dfBarcodes)), collapse = "\n"), "\n\n",
          'File:', fastqFile, '\n',
          'Type:', type, '\n',
          'Output path:', outPath))


  if(type == 'single'){

    outFile <- file.path(outPath, "%.fasta")

    system(paste(fastqMultx,
                 '-x -m 0',
                 '-b', barcodeFile,
                 fastqFile,
                 outFile))
  }

  if(type == 'paired'){

    fileR1 <- fastqFile
    fileR2 <- gsub('_R1', '_R2', fastqFile)
    outR1 <- file.path(outPath, "%_R1.fq.gz")
    outR2 <- file.path(outPath, "%_R2.fq.gz")

    system(paste(fastqMultx,
                 '-x -m 0',
                 '-b', barcodeFile,
                 fileR1, fileR2,
                 "-o", outR1, "-o", outR2))
  }

  unlink(barcodeFile)

  message(paste('Finished demultiplex using:\n',
                'Barcodes:', "\n\n", paste(capture.output(print(dfBarcodes)), collapse = "\n"), "\n\n",
                'File:', fastqFile, '\n',
                'Type:', type, '\n',
                'Output path:', outPath))



}




