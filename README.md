# roplsMini: A lightweight package wrapper of ropls
Removed ExpressionSet, MultiAssayExperiment, MultiDataSet and SummarizedExperiment objects, and only retained PCA, PLS(-DA) and OPLS(-DA) methods for data.frame and matirx in [ropls](https://bioconductor.org/packages/release/bioc/html/ropls.html).

``` r
if (!requireNamespace("remotes", quietly=TRUE))
    install.packages("remotes")
remotes::install_github("SMUZhanLi/roplsMini")

library(roplsMini)

## OPLS-DA
data(sacurine)
dataMatrix <- sacurine$dataMatrix
sampleMetadata <- sacurine$sampleMetadata
variableMetadata <- sacurine$variableMetadata
sacurine.oplsda <- opls(dataMatrix, sampleMetadata[, "gender"], predI = 1, orthoI = NA)
getVipVn(sacurine.oplsda, orthoL = TRUE)
```

