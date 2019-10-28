#' Umi4c preprocessing
#'
#' This function filters the reads keeping the one that present bait sequence, pad and restriction enzyme. 
#' Moreover, it introduces the umi sequence in the header of the fasta files
#' 
#' @inheritParams umi4CatsContacts
#' @export
#' 
#' 
prepScript <- function(input,
                          cores,
                          bait_seq,
                          bait_pad,
                          re,
                          ref_gen,
                          wk_dir,
                          fastqmultx = "fastq-multx"){
# run prep script
  prep_script <- system.file("exec/prep.sh", package = "UMI4Cats")
  system(paste(prep_script,
                "-i", input,
                "-c", cores,
                "-s", bait_seq,
                "-p", bait_pad,
                "-r", re,
                "-g", ref_gen,
                "-w", wk_dir,
                "-f", fastqmultx))
       
}






