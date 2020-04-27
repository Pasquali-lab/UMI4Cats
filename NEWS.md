# UMI4Cats 0.99.2
* Allow `ref_umi4c` to be used as reference for plotting colors, domainogram and differential analysis (not only for normalization).
* Fixed error when using `sampleID` as `grouping` variable in `makeUMI4C()`.
* Fixed bug in `results()` when `fomat=data.frame` and `ordered=TRUE`.

# UMI4Cats 0.99.1
* Fixed error in function `createGeneAnnotation` and `plotGenes` that occurs when there are no genes in the region or a gene has multiple identifiers.
* Fixed duplicated generics definition for `SummarizedExperiment` objects to avoid error when reloading the package.
* Fixed error when `bait_exclusion` is set to 0.
* Added possibility to specify the sample to use as reference for normalization (`ref_umi4c` argument in `makeUMI4C`).
* Now the `grouping` variable in `makeUMI4C()` is used more upstream in the analysis. For using different grouping variables, user must create different `UMI4C` objects.
* Fixed bug where sometimes bait coordinates in the output tsv file are `NA`.
* `statsUMI4C` now also outputs a stats summary table in `wk_dir/logs/stats_summary.txt`.
* Improve function documentation.
* Improve pkgdown UMI4Cats site.
* Rewrite and improve the `Analyzing UMI-4C data with UMI4Cats` vignette.

# UMI4Cats 0.99.0

* Added a `NEWS.md` file to track changes to the package.
* First public release of UMI4Cats.
