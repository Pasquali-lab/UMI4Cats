devtools::load_all()
fastq_dir <- '/imppc/labs/lplab/share/marc/epimutations/raw/fastq/umi4c/prove/fastq_dir'
wk_dir <- '/imppc/labs/lplab/share/marc/epimutations/raw/fastq/umi4c/prove/wk_dir'
bait_seq <- 'ACCTAGAAGGATATGCGCTTGC'
bait_pad <- 'GCGTTAGA'
res_enz <- 'GATC'
cut_seq_5p <- ''
cut_seq_3p <- 'GATC'
ref_gen <- '/imppc/labs/lplab/share/marc/refgen/hg19/hg19.fa'

prepUMI4C(fastq_dir = fastq_dir,
                      wk_dir = wk_dir,
                      bait_seq = bait_seq,
                      bait_pad = bait_pad,
                      res_enz = res_enz)



start_time <- Sys.time()
sleep_for_a_minute()
end_time <- Sys.time()
end_time - start_time
