#' Genomic track generator
#'
#' This function filters the reads keeping the one that present bait sequence, pad and restriction enzyme. 
#' Moreover, it introduces the umi sequence in the header of the fasta files
#' 
#'
#' @param pos Restiction enzyme cut. For example pos 5 in HgaI GACGC indicates cleavage as follows: 5' GACGC^NNNNN.
#' @param re Restriction enzyme sequence.
#' @param nameRe Name restriction enzyme.
#' @param refgen Reference genome used for the generation of the genomic tracks.
#' @param output Output path where to save the genomic tracks.

#' @export


genomicTrackGenerator <- function(pos,
                                  re,
                                  nameRe,
                                  refgen,
                                  output){
  # run prep script
  gtg_script <- system.file("exec/genomicTrackGenerator.sh", package = "UMI4Cats")
  system(paste(gtg_script,
               "-p", pos,
               "-r", re,
               "-n", nameRe,
               "-g", refgen,
               "-o", output))
  
}
