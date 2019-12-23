# UMI4Cats: Processing and analysis of UMI-4C data 

![hexlogo](man/figures/logo.png)

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

## Basic usage

```
library(UMI4Cats)

## 1) Generate Digested genome ----------------------------
# The selected RE in this case is DpnII (|GATC), so the cs5p is "" and cs3p is GATC
hg19_dpnii <- digestGenome(cut_pos = 0,
                           res_enz = "GATC",
                           name_RE = "DpnII",
                           refgen = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, 
                           out_path = "digested_genome/")

## 2) Process UMI-4C fastq files --------------------------
raw_dir <- system.file(file.path("extdata", "SOCS1", "fastq"), 
                       package="UMI4Cats")

contactsUMI4C(fastq_dir = raw_dir,
              wk_dir = "SOCS1",
              bait_seq = "CCCAAATCGCCCAGACCAG",
              bait_pad = "GCGCG",
              res_enz = "GATC",
              cut_pos = 0,
              digested_genome = hg19_dpnii,
              ref_gen = "/biodata/indices/species/Hsapiens/ucsc.hg19.fa",
              threads = 5)
              
## 3) Get filtering and alignment stats -------------------
statsUMI4C(wk_dir = system.file("extdata", "SOCS1",
                               package="UMI4Cats"))

## 4) Analyze UMI-4C results ------------------------------
# Load sample processed file paths
files <- list.files(system.file("extdata", "SOCS1", "count", 
                                package="UMI4Cats"),
                    pattern="*_altCounts.tsv",
                    full.names=TRUE)

# Create colData including all relevant information
colData <- data.frame(sampleID = gsub("_altCounts.tsv", "", basename(files)),
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

## 5) Perform differential test ---------------------------
umi <- fisherUMI4C(umi,
                   filter_low = 20)

## 6) Plot results ----------------------------------------
plotUMI4C(umi, 
          ylim=c(0,10),
          xlim=c(11e6, 11.5e6))
```
