#' Umi4c analysis
#'
#' This function parses, splits, aligns and analyses fastq files counting chromatin contacts between the viewpoint and 
#' the surrounding loci. 
#'
#' @param input Path with the fastq (compressed or descompressed)
#' @param cores Number of cores to use in the analysis
#' @param bait_seq Bait sequence.
#' @param bait_path Bait pad.
#' @param re Restriction enzyme sequence.
#' @param ref_gen Reference genome in fasta format.
#' @param wk_dir Working directory where to carry out the UMI-4c analysis.
#' @param genomic_track Genomic track path
#' @param trimmomatic Trimmomatic path 
#' @param bowtie2 Bowtie2 path 
#' @param fastqmultx Fastq-multx pàth
#' @export


umi4CatsContacts <- function(input,
                             cores,
                             bait_seq,
                             bait_pad,
                             re,
                             ref_gen,
                             wk_dir,
                             genomic_track,
                             trimmomatic="trimmomatic",
                             bowtie2="bowtie2",
                             fastqmultx="fastq-multx"){

  prepScript(input = input,
             cores = cores,
             bait_seq = bait_seq,
             bait_pad = bait_pad,
             re = re,
             ref_gen = ref_gen,
             wk_dir = wk_dir)
  
  splitScript(wk_dir = wk_dir,
              cores = cores,
              re = re,
              trimmomatic = trimmomatic)
  
  alignScript(wk_dir = wk_dir,
              cores = cores,
              ref_gen = ref_gen,
              bowtie2 = bowtie2)
              
              
  umiCounter(wk_dir = wk_dir,
             cores = cores,
             bait_seq = bait_seq,
             bait_pad = bait_pad,
             re = re,
             ref_gen = ref_gen,
             genomic_track = genomic_track,
             bowtie2 = bowtie2)}