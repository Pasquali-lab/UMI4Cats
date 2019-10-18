# UMI4Cats: Processing and analysis of UMI-4C data <img src="man/figures/logo.png" width="121px" height="140px" align="right" style="padding-left:10px;background-color:white;" />

![](https://media.giphy.com/media/JIX9t2j0ZTN9S/giphy-downsized.gif)

## Installation
This is a private repository, so it is a bit tricky to install this package directly 
from R. In order to have a successfull installation you need to follow these steps:

1. Open a bash/shell terminal and clone this repository:
```
git clone https://dogbert.imppc.local/lplab/umi4cats
cd umi4cats
```

2. Open an R/Rstudio session, install the __devtools__ R package and install UMI4Cats.
```
install.packages("devtools")
devtools::install(".")
```

3. Now you can load the package using `library(UMI4Cats)`. 

### Requirements
The R packages requiered for running this package will be automatically checked and installed by devtools. You can check the complete list
in the "Imports" section of the DESCRIPTION file.

However, this package also uses some command line tools that you will need to install manually. Here is the complete list:

- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- [FastqMultx](https://github.com/ExpressionAnalysis/ea-utils/blob/wiki/FastqMultx.md)

## Basic usage

_Work in progress_

![](https://media.giphy.com/media/SU27oG4wJy6Oc/giphy-tumblr.gif)

This pipeline allows the inference of the UMI contacts between a defined viewpoint and the surrounding chromatine given a determinated bait. For each sample-bait combination a UMI-counts table will be generated. 

The main function for the generation of this UMI-counts table is umi4CatsContacts(). This function will analyse all the fastq files saved into a directory with the same bait construction. 

```
input <- '/imppc/labs/lplab/share/marc/umi4cBin/raw/toAnalyse'
cores <- 4
bait_seq <- 'GTTGTCCTTGGGTTTAGCTGC'
bait_pad <- 'ACCTCT'
re <- 'GTAC'
ref_gen <- '/imppc/labs/lplab/msubirana/../share/marc/refgen/hg19/hg19.fa'
wk_dir <- '/imppc/labs/lplab/share/marc/umi4cBin/prove'
genomic_track <- '/imppc/labs/lplab/share/marc/umi4cBin/prove/genomic_tracks_hg19/csp6i_genomicTrack'
trimmomatic <- '/software/debian-8/bio/trimmomatic-0.36/trimmomatic-0.36.jar'
bowtie2 <- 'bowtie2'

umi4CatsContacts(input = input,
                 cores = cores,
                 bait_seq = bait_seq,
                 bait_pad = bait_pad,
                 re = re,
                 ref_gen = ref_gen,
                 wk_dir = wk_dir,
                 genomic_track = genomic_track,
                 trimmomatic = trimmomatic,
                 bowtie2 = bowtie2)
```
