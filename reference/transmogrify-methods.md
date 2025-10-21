# Mogrify a transcriptome using a set of variants

Use a set of SNPs, insertions and deletions to modify a reference
transcriptome

## Usage

``` r
transmogrify(x, var, exons, ...)

# S4 method for class 'XStringSet,GRanges,GRanges'
transmogrify(
  x,
  var,
  exons,
  alt_col = "ALT",
  trans_col = "transcript_id",
  omit_ranges = NULL,
  tag = NULL,
  sep = "_",
  var_tags = FALSE,
  var_sep = "_",
  ol_vars = "fail",
  verbose = TRUE,
  mc.cores = 1,
  ...
)

# S4 method for class 'BSgenome,GRanges,GRanges'
transmogrify(
  x,
  var,
  exons,
  alt_col = "ALT",
  trans_col = "transcript_id",
  omit_ranges = NULL,
  tag = NULL,
  sep = "_",
  var_tags = FALSE,
  var_sep = "_",
  ol_vars = "fail",
  verbose = TRUE,
  mc.cores = 1,
  ...
)

# S4 method for class 'BSgenome,VcfFile,GRanges'
transmogrify(
  x,
  var,
  exons,
  alt_col = "ALT",
  trans_col = "transcript_id",
  omit_ranges = NULL,
  tag = NULL,
  sep = "_",
  var_tags = FALSE,
  var_sep = "_",
  ol_vars = "fail",
  verbose = TRUE,
  mc.cores = 1,
  which,
  ...
)

# S4 method for class 'XStringSet,VcfFile,GRanges'
transmogrify(
  x,
  var,
  exons,
  alt_col = "ALT",
  trans_col = "transcript_id",
  omit_ranges = NULL,
  tag = NULL,
  sep = "_",
  var_tags = FALSE,
  var_sep = "_",
  ol_vars = "fail",
  verbose = TRUE,
  mc.cores = 1,
  which,
  ...
)
```

## Arguments

- x:

  Reference genome as either a DNAStringSet or BSgenome

- var:

  GRanges object containing the variants

- exons:

  GRanges object with ranges representing exons

- ...:

  Passed to
  [parallel::mclapply](https://rdrr.io/r/parallel/mclapply.html)

- alt_col:

  Column from `var` containing alternate bases

- trans_col:

  Column from 'exons' containing the transcript_id

- omit_ranges:

  GRanges object containing ranges to omit, such as PAR-Y regions, for
  example

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

  logical(1) Include informative messages, or operate silently

- mc.cores:

  Number of cores to be used when multi-threading via
  [parallel::mclapply](https://rdrr.io/r/parallel/mclapply.html)

- which:

  GRanges object passed to
  [VariantAnnotation::ScanVcfParam](https://rdrr.io/pkg/VariantAnnotation/man/ScanVcfParam-class.html)
  if using a VCF directly

## Value

An XStringSet

## Details

Produce a set of variant modified transcript sequences from a standard
reference genome. Supported variants are SNPs, Insertions and Deletions

Ranges needing to be masked, such as the Y-chromosome, or Y-PAR can be
provided.

**It should be noted that this is a time consuming process** Inclusion
of a large set of insertions and deletions across an entire
transcriptome can involve individually modifying many thousands of
transcripts, which can be a computationally demanding task. Whilst this
can be parallelised using an appropriate number of cores, this may also
prove taxing for lower power laptops, and pre-emptively closing memory
hungry programs such as Slack, or internet browers may be prudent.

## Examples

``` r
library(GenomicRanges)
library(GenomicFeatures)
seq <- DNAStringSet(c(chr1 = "ACGTAAATGG"))
exons <- GRanges(c("chr1:1-3:-", "chr1:7-9:-"))
exons$transcript_id <- c("trans1")

# When using extractTranscriptSeqs -stranded exons need to be sorted by end
exons <- sort(exons, decreasing = TRUE, by = ~end)
exons
#> GRanges object with 2 ranges and 1 metadata column:
#>       seqnames    ranges strand | transcript_id
#>          <Rle> <IRanges>  <Rle> |   <character>
#>   [1]     chr1       7-9      - |        trans1
#>   [2]     chr1       1-3      - |        trans1
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
trByExon <- splitAsList(exons, exons$transcript_id)

# Check the sequences
seq
#> DNAStringSet object of length 1:
#>     width seq                                               names               
#> [1]    10 ACGTAAATGG                                        chr1
extractTranscriptSeqs(seq, trByExon)
#> DNAStringSet object of length 1:
#>     width seq                                               names               
#> [1]     6 CATCGT                                            trans1

# Define some variants
var <- GRanges(c("chr1:2", "chr1:8"))
var$ALT <- c("A", "GGG")

# Include the variants adding tags to indicate a SNP and indel
# The exons GRanges object will be split by transcript internally
transmogrify(seq, var, exons, var_tags = TRUE)
#> 1 transcripts found with indels
#> DNAStringSet object of length 1:
#>     width seq                                               names               
#> [1]     8 CCCCTCTT                                          trans1_si

```
