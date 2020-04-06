# UMI4Cats 0.99.1
* Fixed error in function `createGeneAnnotation` and `plotGenes` that occurs when there are no genes in the region or a gene has multiple identifiers.
* Fixed duplicated generics definition for `SummarizedExperiment` objects to avoid error when reloading the package.
* Fixed error when `bait_exclusion` is set to 0.
* Added possibility to specify the sample to use as reference for normalization (`ref_umi4c` argument in `makeUMI4C`).
* Now the `grouping` variable in `makeUMI4C()` is used more upstream in the analysis. For using different grouping variables, user must create different `UMI4C` objects.
* Fixed bug where sometimes bait coordinates in the output tsv file are `NA`.
* Minor changes in the package vignette. 

# UMI4Cats 0.99.0

* Added a `NEWS.md` file to track changes to the package.
* First public release of UMI4Cats.
