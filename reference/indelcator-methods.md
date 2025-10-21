# Substitute InDels into one or more sequences

Modify one or more sequences to include Insertions or Deletions

## Usage

``` r
indelcator(x, indels, ...)

# S4 method for class 'XString,GRanges'
indelcator(x, indels, exons, alt_col = "ALT", ol_vars = "fail", ...)

# S4 method for class 'DNAStringSet,GRanges'
indelcator(
  x,
  indels,
  alt_col = "ALT",
  ol_vars = "fail",
  mc.cores = 1,
  verbose = TRUE,
  ...
)

# S4 method for class 'BSgenome,GRanges'
indelcator(
  x,
  indels,
  alt_col = "ALT",
  ol_vars = "fail",
  mc.cores = 1,
  names,
  ...
)
```

## Arguments

- x:

  Sequence of class XString

- indels:

  GRanges object with InDel locations and the alternate allele

- ...:

  Passed to
  [parallel::mclapply](https://rdrr.io/r/parallel/mclapply.html)

- exons:

  GRanges object containing exon structure for `x`

- alt_col:

  Column containing the alternate allele

- ol_vars:

  Error handling for any overlapping variants. See
  [cleanVariants](https://smped.github.io/transmogR/reference/cleanVariants-methods.md)
  for possible values and an explanation

- mc.cores:

  Number of cores to use when calling
  [parallel::mclapply](https://rdrr.io/r/parallel/mclapply.html)
  internally

- verbose:

  logical(1) Print all messages

- names:

  passed to
  [`BSgenome::getSeq()`](https://rdrr.io/pkg/Biostrings/man/getSeq.html)
  when x is a BSgenome object

## Value

A DNAStringSet or XString object (See Details)

## Details

This is a lower-level function relied on by both
[`transmogrify()`](https://smped.github.io/transmogR/reference/transmogrify-methods.md)
and
[`genomogrify()`](https://smped.github.io/transmogR/reference/genomogrify-methods.md).

Takes an
[Biostrings::XString](https://rdrr.io/pkg/Biostrings/man/XString-class.html)
or
[Biostrings::XStringSet](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
object and modifies the sequence to incorporate InDels. The expected
types of data determine the behaviour, with the following expectations
describing how the function will incorporate data

|  |  |  |  |
|----|----|----|----|
| Input Data Type | Exons Required | Use Case | Returned |
| XString | Y | Modify a Reference Transcriptome | XString |
| DNAStringSet | N | Modify a Reference Genome | DNAStringSet |
| BSgenome | N | Modify a Reference Genome | DNAStringSet |

## See also

[`transmogrify()`](https://smped.github.io/transmogR/reference/transmogrify-methods.md)
[`genomogrify()`](https://smped.github.io/transmogR/reference/genomogrify-methods.md)

## Examples

``` r
## Start with a DNAStringSet
library(GenomicRanges)
seq <- DNAStringSet(c(seq1 = "AATCTGCGC"))
## Define an Insertion
var <- GRanges("seq1:1")
var$ALT <- "AAA"
seq
#> DNAStringSet object of length 1:
#>     width seq                                               names               
#> [1]     9 AATCTGCGC                                         seq1
indelcator(seq, var)
#> Updating seq1; Original length: 9
#> ; Updated length: 11
#> DNAStringSet object of length 1:
#>     width seq                                               names               
#> [1]    11 AAAATCTGCGC                                       seq1

## To modify a single transcript
library(GenomicFeatures)
#> Loading required package: AnnotationDbi
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
ex <- GRanges(c("seq1:1-3:+", "seq1:7-9:+"))
orig <- extractTranscriptSeqs(seq, GRangesList(tx1 = ex))[["tx1"]]
orig
#> 6-letter DNAString object
#> seq: AATCGC
indelcator(orig, var, exons = ex)
#> 8-letter DNAString object
#> seq: AAAATCGC
```
