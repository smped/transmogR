# Get the PAR-Y Regions From a Seqinfo Object

Define the Pseudo-Autosomal Regions from a Seqinfo Object

## Usage

``` r
parY(x, ...)

# S4 method for class 'Seqinfo'
parY(x, ...)

# S4 method for class 'character'
parY(x, prefix = NULL, ...)
```

## Arguments

- x:

  A Seqinfo object or any of named build. If passing a character vector,
  [`match.arg()`](https://rdrr.io/r/base/match.arg.html) will be used to
  match the build.

- ...:

  Not used

- prefix:

  Optional prefix to place before chromosome names. Can only be NULL, ""
  or "chr"

## Value

A GenomicRanges object

## Details

Using a seqinfo object based on either hg38, hg19, CHM13.v2 or their
variations, create a GRanges object with the Pseudo-Autosomal Regions
from the Y chromosome for that build. The length of the Y chromosome on
the seqinfo object is used to determine the correct genome build when
passing a Seqinfo object. Otherwise

An additional mcols column called PAR will indicate PAR1 and PAR2

## Examples

``` r
library(Seqinfo)
sq <- Seqinfo(
    seqnames = "chrY", seqlengths = 59373566, genome = "hg19_only_chrY"
)
parY(sq)
#> GRanges object with 2 ranges and 1 metadata column:
#>       seqnames            ranges strand |         PAR
#>          <Rle>         <IRanges>  <Rle> | <character>
#>   [1]     chrY     10001-2649520      * |        PAR1
#>   [2]     chrY 59034050-59363566      * |        PAR2
#>   -------
#>   seqinfo: 1 sequence from hg19_only_chrY genome

## PAR regions for CHM13 are also available
sq <- Seqinfo(
    seqnames = "chrY", seqlengths = 62460029, genome = "CHM13"
)
parY(sq)
#> GRanges object with 2 ranges and 1 metadata column:
#>       seqnames            ranges strand |         PAR
#>          <Rle>         <IRanges>  <Rle> | <character>
#>   [1]     chrY         1-2458320      * |        PAR1
#>   [2]     chrY 62122809-62460029      * |        PAR2
#>   -------
#>   seqinfo: 1 sequence from CHM13 genome

## Or just call by name
parY("GRCh38", prefix = "chr")
#> GRanges object with 2 ranges and 1 metadata column:
#>       seqnames            ranges strand |         PAR
#>          <Rle>         <IRanges>  <Rle> | <character>
#>   [1]     chrY     10001-2781479      * |        PAR1
#>   [2]     chrY 56887903-57217415      * |        PAR2
#>   -------
#>   seqinfo: 1 sequence from GRCh38 genome; no seqlengths

```
