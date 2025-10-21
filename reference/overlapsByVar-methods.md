# Count overlaps by variant type

Count how many variants of each type overlap ranges

## Usage

``` r
overlapsByVar(x, var, ...)

# S4 method for class 'GRangesList,GRanges'
overlapsByVar(x, var, alt_col = "ALT", ...)

# S4 method for class 'GRanges,GRanges'
overlapsByVar(x, var, alt_col = "ALT", ...)
```

## Arguments

- x:

  A GRangesList with features of interest

- var:

  A Granges object with variants of interest

- ...:

  Passed to
  [rowSums](https://rdrr.io/pkg/MatrixGenerics/man/rowSums.html)

- alt_col:

  The column within mcols(var) which contains the alternate allele

## Value

A vector or matrix

## Details

Taking any GRanges or GRangesList, count how many of each variant type
overlap a region.

## Examples

``` r
library(rtracklayer)
library(VariantAnnotation)
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: ‘matrixStats’
#> The following objects are masked from ‘package:Biobase’:
#> 
#>     anyMissing, rowMedians
#> 
#> Attaching package: ‘MatrixGenerics’
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> The following object is masked from ‘package:Biobase’:
#> 
#>     rowMedians
#> Loading required package: SummarizedExperiment
#> Loading required package: Rsamtools
#> 
#> Attaching package: ‘VariantAnnotation’
#> The following object is masked from ‘package:base’:
#> 
#>     tabulate
gtf <- import.gff(
    system.file("extdata/gencode.v44.subset.gtf.gz", package = "transmogR")
)
grl <- splitAsList(gtf, gtf$type)
vcf <- system.file("extdata/1000GP_subset.vcf.gz", package = "transmogR")
var <- rowRanges(readVcf(vcf, param = ScanVcfParam(fixed = "ALT")))
overlapsByVar(grl, var)
#>            Deletion Insertion SNV Total
#> gene              6         1  91    98
#> transcript        6         1  91    98
#> exon              0         1  37    38
```
