# Create a set of tags indicating overlap status with variants

Create a set of tags indicating overlap status with variants

## Usage

``` r
varTags(x, var, tag = NULL, var_tags = TRUE, sep = "_", pre = sep, ...)
```

## Arguments

- x:

  GRanges or GRangesList

- var:

  Set of variants for `x` to be compared to

- tag:

  Tag to be added for all overlapping positions

- var_tags:

  logical(1) Include 's', 'i' and 'd' tags. See details

- sep:

  Separator added between tag and var_tags

- pre:

  Separator to add at the start of returned tags

- ...:

  Passed to
  [`cleanVariants()`](https://smped.github.io/transmogR/reference/cleanVariants-methods.md)

## Value

Character vector of the same length as `x`

## Details

Take a GRanges or GRangesList and compare against a set of variants.
Variants will be classified into SNV, Insertions and Deletions using
[`varTypes()`](https://smped.github.io/transmogR/reference/varTypes.md)
and tags defined. An overall set of tags defining any overlap can be
created by themselves. An additional set of tags containing 's', 'i' or
'd' to indicate overlap with an SNV, Insertion or Deletion can also be
created, with the concatentation of both tags being returned.

## Examples

``` r
# Load the included subset of 1000 Genomes Variants
library(VariantAnnotation)
vcf <- system.file("extdata/1000GP_subset.vcf.gz", package = "transmogR")
vcf <- VcfFile(vcf)
var <- cleanVariants(vcf)
# Now load some exons, then split by transcript, subsetting to the first 40
library(rtracklayer)
f <- system.file("extdata/gencode.v44.subset.gtf.gz", package = "transmogR")
gtf <- import.gff(f, feature.type = "exon")
exon_by_trans <- splitAsList(gtf, gtf$transcript_id)[1:40]
# And produce tags based on the overlapping variants within the exons
# Overlapping SNVs will return an 's' whilst insertions include an 'i'
varTags(exon_by_trans, var, tag = "1000GP")
#>  [1] "_1000GP_s"  ""           "_1000GP_s"  ""           ""          
#>  [6] ""           ""           ""           "_1000GP_s"  ""          
#> [11] ""           "_1000GP_s"  "_1000GP_s"  ""           ""          
#> [16] "_1000GP_s"  ""           ""           ""           ""          
#> [21] "_1000GP_s"  ""           "_1000GP_s"  ""           ""          
#> [26] ""           ""           ""           ""           ""          
#> [31] "_1000GP_s"  ""           ""           "_1000GP_si" ""          
#> [36] "_1000GP_s"  ""           ""           ""           "_1000GP_s" 

```
