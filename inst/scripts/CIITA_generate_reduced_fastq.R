######################################
# Generate reduced CIITA fastq files #
######################################
# This script describes how files contained in the following folders were 
# generated.
#    - inst/extdata/CIITA/fastq
#
# Briefly, CIITA raw fastq files were downloaded from GEO, demultiplexed
# according to the samples of origin and subsampled to reduce file size and
# thus computation time for downstream analyses.

# 1) Download raw data from source: --------------------------------------------
#   - CIITA control: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4059911
#   - CIITA cytokines: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4059926

# 2) Demultiplex different samples: --------------------------------------------
files <- list.files("~/Downloads/",
                    pattern = ".fastq.gz.1",
                    full.names = TRUE)

# Prepare output files
out_files_hi24 <- gsub("CIITA_m", "hi24_CIITA", 
                       gsub("fastq.gz.1", "fastq.gz", 
                            basename(files)))

out_files_hi32 <- gsub("CIITA_m", "hi32_CIITA", 
                       gsub("fastq.gz.1", "fastq.gz", 
                            basename(files)))

# Prepare filtering functions
filterFun_hi24 <- function(reads) {
  reads[as.character(subseq(id(reads), start=22, end=22)) == 5]
}

filterFun_hi32 <- function(reads) {
  reads[as.character(subseq(id(reads), start=22, end=22)) == 6]
}

# Demultiplex files
ShortRead::filterFastq(files,
                       destinations = out_files_hi24,
                       filter = filterFun_hi24)

ShortRead::filterFastq(files,
                       destinations = out_files_hi32,
                       filter = filterFun_hi32)

# 3) Subsample fastq files ----------------------------------------------------
# Fastq files in inst/extdata: 100 total reads
# Fastq files in downloadUMI4CexampleData: 200K total reads
# Change following argument to sample the desired number of reads

# reads <- 2e5 # downloadUMI4CexampleData
reads <- 100 # inst/extdata

fastq_files <- c(out_files_hi24, 
                 out_files_hi32)

out_dir <- paste0("sub_", reads)
dir.create(out_dir, FALSE)

for (file in fastq_files) {
  message(">> ", file)
  
  out_file <- file.path(out_dir, file)

  cmd <- paste("seqtk sample",
               "-s123",
               file,
               reads,
               "| gzip -c",
               ">", out_file)
  system(cmd)
  print(cmd)
}