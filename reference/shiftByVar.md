# Calculate new exon co-ordinates

Calcluate new exon co-ordinates after including InDels

## Usage

``` r
shiftByVar(x, var, alt_col = "ALT", mc.cores = 1, ...)
```

## Arguments

- x:

  A GenomicRanges object with co-ordinates needing to be recalculated

- var:

  A set of variants to be incorporated into a reference genome

- alt_col:

  The name of the column with the alternate sequence for each variant

- mc.cores:

  Passed internally to
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)

- ...:

  Not used

## Value

GRanges object with co-ordinates shifted according the to the provided
variants. The new co-ordinates will be compatible with a
variant-modified genome as produced by
[`genomogrify()`](https://smped.github.io/transmogR/reference/genomogrify-methods.md)
and can be used to extract the sequences associated with the ranges in
the modified reference.

## Details

Given a set of variants, this will return a set of genomic ranges with
updated co-ordinates able to be applied on a variant modified reference
genome

## Examples

``` r
# Define a 3nt insertion
var <- GRanges("seq1:5:*", seqlengths = c(seq1=10), REF = "A", ALT = "AGT")
var
#> GRanges object with 1 range and 2 metadata columns:
#>       seqnames    ranges strand |         REF         ALT
#>          <Rle> <IRanges>  <Rle> | <character> <character>
#>   [1]     seq1         5      * |           A         AGT
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
# A simple GRanges to shift co-ordinates for
gr <- GRanges("seq1:1-10:+", seqlengths = c(seq1=10), feature = "feature1")
gr
#> GRanges object with 1 range and 1 metadata column:
#>       seqnames    ranges strand |     feature
#>          <Rle> <IRanges>  <Rle> | <character>
#>   [1]     seq1      1-10      + |    feature1
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
# Create shifted co-ordinates based on the provided variants
new_gr <- shiftByVar(gr, var)
new_gr
#> GRanges object with 1 range and 1 metadata column:
#>       seqnames    ranges strand |     feature
#>          <Rle> <IRanges>  <Rle> | <character>
#>   [1]     seq1      1-12      + |    feature1
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
## The seqlengths will have been adjusted to account for all variants
seqinfo(new_gr)
#> Seqinfo object with 1 sequence from an unspecified genome:
#>   seqnames seqlengths isCircular genome
#>   seq1             12         NA   <NA>
```
