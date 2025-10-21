# OverWrite Letters in an XStringSet

OverWrite Letters (e.g. SNPs) in an XStringSet

## Usage

``` r
owl(seq, snps, ...)

# S4 method for class 'XStringSet,GRanges'
owl(seq, snps, alt_col = "ALT", ol_vars = "fail", ...)

# S4 method for class 'BSgenome,GRanges'
owl(seq, snps, alt_col = "ALT", names, ...)
```

## Arguments

- seq:

  A BSgenome, DNAStringSet, RNAStringSet or other XStringSet.

- snps:

  A GRanges object with SNP positions and a column containing the
  alternate allele

- ...:

  Passed to
  [`Biostrings::replaceLetterAt()`](https://rdrr.io/pkg/Biostrings/man/replaceLetterAt.html)

- alt_col:

  Column name in the mcols element of `snps` containing the alternate
  allele

- ol_vars:

  Error handling for any overlapping variants. See
  [cleanVariants](https://smped.github.io/transmogR/reference/cleanVariants-methods.md)
  for possible values and an explanation

- names:

  Sequence names to operate on

## Value

An object of the same class as the original object, but with SNPs
inserted at the supplied positions

## Details

This is a lower-level function called by
[`transmogrify()`](https://smped.github.io/transmogR/reference/transmogrify-methods.md)
and
[`genomogrify()`](https://smped.github.io/transmogR/reference/genomogrify-methods.md),
but able to be called by the user if needed

Note that when providing a BSgenome object, this will first be coerced
to a DNAStringSet which can be time consuming.

## Examples

``` r
seq <- DNAStringSet(c(chr1 = "AAGC"))
snps <- GRanges("chr1:2")
snps$ALT <- "G"
snps
#> GRanges object with 1 range and 1 metadata column:
#>       seqnames    ranges strand |         ALT
#>          <Rle> <IRanges>  <Rle> | <character>
#>   [1]     chr1         2      * |           G
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
seq
#> DNAStringSet object of length 1:
#>     width seq                                               names               
#> [1]     4 AAGC                                              chr1
owl(seq, snps)
#> DNAStringSet object of length 1:
#>     width seq                                               names               
#> [1]     4 AGGC                                              chr1
```
