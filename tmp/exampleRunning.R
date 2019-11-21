fastq_dir <- '/imppc/labs/lplab/share/marc/epimutations/raw/fastq/umi4c/prove/fastq_dir'
wk_dir <- '/imppc/labs/lplab/share/marc/epimutations/raw/fastq/umi4c/prove/wk_dir'
bait_seq <- 'ACCTAGAAGGATATGCGCTTGC'
bait_pad <- 'GCGTTAGA'
res_enz <- 'GATC'

prepUMI4C(fastq_dir = fastq_dir,
                      wk_dir = wk_dir,
                      bait_seq = bait_seq,
                      bait_pad = bait_pad,
                      res_enz = res_enz)
