# UMI4Cats: Processing and analysis of UMI-4C data

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
- [Picard Tools](https://broadinstitute.github.io/picard/)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

## Basic usage

_Work in progress_

![](https://media.giphy.com/media/SU27oG4wJy6Oc/giphy-tumblr.gif)
