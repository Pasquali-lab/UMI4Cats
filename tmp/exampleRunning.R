devtools::load_all()
fastq_dir <- '/imppc/labs/lplab/share/marc/bCell/raw/fastq/umi4c/csp6i_CCR1-5_TFRC_enh/S1'
wk_dir <- '/imppc/labs/lplab/share/marc/bCell/processed/hg19/umi4c/r'
bait_seq <- 'GTTGTCCTTGGGTTTAGCTGC'
bait_pad <- 'ACCTCT'
res_enz <- 'GTAC'
cut_pos <- 1
ref_gen <- '/imppc/labs/lplab/share/marc/refgen/hg19/hg19.fa'
threads <- 6
out_path <- '/imppc/labs/lplab/share/marc/refgen/digestGenome/'
name_RE <- 'csp6i'
library(BSgenome.Hsapiens.UCSC.hg19)
ref_gen <- BSgenome.Hsapiens.UCSC.hg19
digested_genome <- '/imppc/labs/lplab/share/marc/refgen/digestGenome/BSgenome.Hsapiens.UCSC.hg19_csp6i.tsv'

prepUMI4C(fastq_dir = fastq_dir,
          wk_dir = wk_dir,
          bait_seq = bait_seq,
          bait_pad = bait_pad,
          res_enz = res_enz)

splitUMI4C(wk_dir = wk_dir,
           prep_dir = prep_dir,
           res_enz = res_enz,
           cut_pos = cut_pos)

alignmentUMI4C(wk_dir = wk_dir,
               bait_seq = bait_seq,
               bait_pad = bait_pad,
               res_enz = res_enz,
               ref_gen = ref_gen,
               threads = threads)

counterUMI4C(wk_dir=wk_dir,
             bait_seq=bait_seq,
             bait_pad=bait_pad,
             res_enz=res_enz,
             digested_genome=digested_genome,
             ref_gen=ref_gen)

digestGenome(res_enz = res_enz,
             cut_pos = cut_pos,
             name_RE = name_RE,
             ref_gen = ref_gen,
             out_path = out_path)




########################



# Load sample processed file paths
path <- '/imppc/labs/lplab/share/marc/bCell/processed/hg19/umi4c/r/count'

files <- list.files(path,
                    pattern="*.tsv",
                    full.names=TRUE)

# Create colData including all relevant information
colData <- data.frame(sampleID = gsub("_counts.tsv", "", basename(files)),
                      file = files,
                      stringsAsFactors=F)

library(tidyr)
colData <- colData %>%
  separate(sampleID,
           into=c("condition", "replicate", "viewpoint"),
           remove=FALSE)

# Load UMI-4C data and generate UMI4C object
umi <- makeUMI4C(colData=colData,
                 viewpoint_name="SOCS1")

# Perform differential test
umi <- fisherUMI4C(umi)

# Plot results
plotUMI4C(umi, ylim=c(0,4))
