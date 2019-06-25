#' Umi4c alignment
#'
#' This function aligns the splited reads using a defined reference genome.  
#'
#' @inheritParams umi4CatsContacts
#' @export
#' 
alignScript <- function(wk_dir,
                        cores,
                        ref_gen,
                        bowtie2="bowtie2"){
  # run prep script
  align_script <- system.file("exec/alignment.sh", package = "UMI4Cats")
  system(paste(align_script,
               "-w", wk_dir,
               "-c", cores,
               "-i", ref_gen,
               "-b", bowtie2))
  
}