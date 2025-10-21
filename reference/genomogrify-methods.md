# Mogrify a genome using a set of variants

Use a set of SNPS, insertions and deletions to modify a reference genome

## Usage

``` r
genomogrify(x, var, ...)

# S4 method for class 'XStringSet,GRanges'
genomogrify(
  x,
  var,
  alt_col = "ALT",
  mask = GRanges(),
  tag = NULL,
  sep = "_",
  var_tags = FALSE,
  var_sep = "_",
  ol_vars = "fail",
  verbose = TRUE,
  ...
)

# S4 method for class 'BSgenome,GRanges'
genomogrify(
  x,
  var,
  alt_col = "ALT",
  mask = GRanges(),
  names,
  tag = NULL,
  sep = "_",
  var_tags = FALSE,
  var_sep = "_",
  ol_vars = "fail",
  verbose = TRUE,
  ...
)

# S4 method for class 'BSgenome,VcfFile'
genomogrify(
  x,
  var,
  alt_col = "ALT",
  mask = GRanges(),
  names,
  tag = NULL,
  sep = "_",
  var_tags = FALSE,
  var_sep = "_",
  ol_vars = "fail",
  which,
  verbose = TRUE,
  ...
)

# S4 method for class 'XStringSet,VcfFile'
genomogrify(
  x,
  var,
  alt_col = "ALT",
  mask = GRanges(),
  tag = NULL,
  sep = "_",
  var_tags = FALSE,
  var_sep = "_",
  ol_vars = "fail",
  which,
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  A DNAStringSet or BSgenome

- var:

  GRanges object containing the variants, or a
  [VariantAnnotation::VcfFile](https://rdrr.io/pkg/VariantAnnotation/man/VcfFile-class.html)

- ...:

  Passed to
  [parallel::mclapply](https://rdrr.io/r/parallel/mclapply.html)

- alt_col:

  The name of the column with `var` containing alternate bases

- mask:

  Optional GRanges object defining regions to be masked with an 'N'

- tag:

  Optional tag to add to all sequence names which were modified

- sep:

  Separator to place between seqnames names & tag

- var_tags:

  logical(1) Add tags indicating which type of variant were
  incorporated, with 's', 'i' and 'd' representing SNPs, Insertions and
  Deletions respectively

- var_sep:

  Separator between any previous tags and variant tags

- ol_vars:

  Error handling for any overlapping variants. See
  [cleanVariants](https://smped.github.io/transmogR/reference/cleanVariants-methods.md)
  for possible values and an explanation

- verbose:

  logical(1) Print progress messages while running

- names:

  Sequence names to be mogrified

- which:

  GRanges object passed to
  [VariantAnnotation::ScanVcfParam](https://rdrr.io/pkg/VariantAnnotation/man/ScanVcfParam-class.html)
  if using a VCF directly

## Value

XStringSet with variant modified sequences

## Details

This function is designed to create a variant-modified reference genome,
intended to be included as a set of decoys when using salmon in
selective alignment mode. Sequence lengths will change if InDels are
included and any coordinate-based information will be lost on the output
of this function.

Tags are able to be added to any modified sequence to assist identifying
any changes that have been made to a sequence.

## Examples

``` r
library(GenomicRanges)
dna <- DNAStringSet(c(chr1 = "ACGT", chr2 = "AATTT"))
var <- GRanges(c("chr1:1", "chr1:3", "chr2:1-3"))
var$ALT <- c("C", "GG", "A")
dna
#> DNAStringSet object of length 2:
#>     width seq                                               names               
#> [1]     4 ACGT                                              chr1
#> [2]     5 AATTT                                             chr2
genomogrify(dna, var)
#> Updating chr1; Original length: 4
#> ; Updated length: 5
#> Updating chr2; Original length: 5
#> ; Updated length: 3
#> DNAStringSet object of length 2:
#>     width seq                                               names               
#> [1]     5 CCGGT                                             chr1
#> [2]     3 ATT                                               chr2
genomogrify(dna, var, tag = "mod")
#> Updating chr1; Original length: 4
#> ; Updated length: 5
#> Updating chr2; Original length: 5
#> ; Updated length: 3
#> DNAStringSet object of length 2:
#>     width seq                                               names               
#> [1]     5 CCGGT                                             chr1_mod
#> [2]     3 ATT                                               chr2_mod
genomogrify(dna, var, var_tags = TRUE)
#> Updating chr1; Original length: 4
#> ; Updated length: 5
#> Updating chr2; Original length: 5
#> ; Updated length: 3
#> DNAStringSet object of length 2:
#>     width seq                                               names               
#> [1]     5 CCGGT                                             chr1_si
#> [2]     3 ATT                                               chr2_d
genomogrify(dna, var, mask = GRanges("chr2:1-5"), var_tags = TRUE)
#> Applying mask
#> Updating chr1; Original length: 4
#> ; Updated length: 5
#> DNAStringSet object of length 2:
#>     width seq                                               names               
#> [1]     5 CCGGT                                             chr1_si
#> [2]     5 NNNNN                                             chr2

```
