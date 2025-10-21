# Check provided variants for compatibility

Check and cleanup a set of variants for downstream compatibility

## Usage

``` r
cleanVariants(var, ...)

# S4 method for class 'GRanges'
cleanVariants(var, ol_vars = "fail", ref_col = "REF", alt_col = "ALT", ...)

# S4 method for class 'VcfFile'
cleanVariants(
  var,
  ol_vars = "fail",
  which,
  ref_col = "REF",
  alt_col = "ALT",
  ...
)
```

## Arguments

- var:

  GRanges object containing the variants, or a
  [VariantAnnotation::VcfFile](https://rdrr.io/pkg/VariantAnnotation/man/VcfFile-class.html)

- ...:

  Not used

- ol_vars:

  Error handling for any overlapping variants. Can take values in
  c("fail", "none", "first", "last", "longest", "shortest"). Default is
  set to fail, with additional options to drop all overlapping variants
  ('none'), select by genomic position ('first', 'last'), or select by
  the scale of change to the genome ('longest', 'shortest')

- ref_col, alt_col:

  Column names corresponding to the reference and alternate alleles

- which:

  Passed to
  [VariantAnnotation::ScanVcfParam](https://rdrr.io/pkg/VariantAnnotation/man/ScanVcfParam-class.html)
  if working with a VcfFile

## Value

GRanges object with any incompatible variants removed, or an error
produced. The mcols will contain the columns REF and ALT, unless
otherwise specified, as character vectors

## Details

This function checks a set of variants for the expected structure which
is required by all downstream functions in transmogR. The primary change
to the data structure is that both REF and ALT columns will be set as
simple character vectors.

Given the complicated variant calls that can often be produced by
variant callers, additional checks performed will be to ensure that:

- there are no overlapping variants

- SNPs are all single nucleotides and not longer strings

- SNPs are bi-allelic

- Insertions contain a single nucleotide in the REF column

- Deletions contain a single nucleotide in the ALT column

- No missing values are present

- All ALT/REF nucleotides conform to IUPAC codings

## Examples

``` r
# Any conflicting variants will be removed
var <- GRanges(c("chr10:114468420-114468422", "chr10:114468422"))
var$REF <- c("GCC", "C")
var$ALT <- c("G", "CTAT")
var
#> GRanges object with 2 ranges and 2 metadata columns:
#>       seqnames              ranges strand |         REF         ALT
#>          <Rle>           <IRanges>  <Rle> | <character> <character>
#>   [1]    chr10 114468420-114468422      * |         GCC           G
#>   [2]    chr10           114468422      * |           C        CTAT
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
## Taken from the 1000GP, the first variant would delete the C at 114468422
## whilst the second variant begins an insertion at this position.
## These are clearly conflicting. The default value for ol_vars is to fail
## with an error (ol_vars = "fail"). However, both can be removed by setting
## ol_vars = "none". A warning will always be produced.
cleanVariants(var, ol_vars = "none")
#> Warning: 2 pairs of overlapping loci found and cannot be incorporated into a modified reference.
#> All overlapping variants will be removed
#> GRanges object with 0 ranges and 2 metadata columns:
#>    seqnames    ranges strand |         REF         ALT
#>       <Rle> <IRanges>  <Rle> | <character> <character>
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
## Or the longest can be retained, along with multiple other options
cleanVariants(var, ol_vars = "longest")
#> Warning: 2 pairs of overlapping loci found and cannot be incorporated into a modified reference.
#> The longest overlapping locus by genomic change will be retained
#> GRanges object with 1 range and 2 metadata columns:
#>       seqnames    ranges strand |         REF         ALT
#>          <Rle> <IRanges>  <Rle> | <character> <character>
#>   [1]    chr10 114468422      * |           C        CTAT
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
