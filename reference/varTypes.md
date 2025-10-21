# Identify SNVs, Insertions and Deletions

Identify SNVs, Insertions and Deletions within a GRanges object

## Usage

``` r
varTypes(x, alt_col = "ALT", ...)
```

## Arguments

- x:

  GenomicRanges object

- alt_col:

  Name of the column with mcols(x) which contains the alternate allele.
  Can be an XStringSetList, XStringSet or character

- ...:

  Not used

## Value

Character vector

## Details

Using the width of the reference and alternate alleles, classify each
range as an SNV, Insertion or Deletion.

- SNVs are expected to have REF & ALT widths of 1

- Insertions are expected to have ALT longer than REF

- Deletions are expected to have ALT shorter than REF

These are relatively permissive criteria

## Examples

``` r
# Load the example VCF and classify ranges
library(VariantAnnotation)
f <- system.file("extdata/1000GP_subset.vcf.gz", package = "transmogR")
vcf <- readVcf(f)
gr <- rowRanges(vcf)
type <- varTypes(gr)
table(type)
#> type
#>  Deletion Insertion       SNV 
#>         6         1        93 
gr[type != "SNV"]
#> GRanges object with 7 ranges and 5 metadata columns:
#>                         seqnames        ranges strand | paramRangeID
#>                            <Rle>     <IRanges>  <Rle> |     <factor>
#>          1:788418:CAG:C     chr1 788418-788420      * |           NA
#>       1:789568:TATGGA:T     chr1 789568-789573      * |           NA
#>   1:790933:CGAATGGAAT:C     chr1 790933-790943      * |           NA
#>          1:818464:CCT:C     chr1 818464-818466      * |           NA
#>           1:826577:A:AT     chr1        826577      * |           NA
#>   1:839480:CACACACCTG:C     chr1 839480-839494      * |           NA
#>   1:839515:CTAGACACAC:C     chr1 839515-839543      * |           NA
#>                                             REF                ALT      QUAL
#>                                  <DNAStringSet> <DNAStringSetList> <numeric>
#>          1:788418:CAG:C                     CAG                  C        NA
#>       1:789568:TATGGA:T                  TATGGA                  T        NA
#>   1:790933:CGAATGGAAT:C             CGAATGGAATG                  C        NA
#>          1:818464:CCT:C                     CCT                  C        NA
#>           1:826577:A:AT                       A                 AT        NA
#>   1:839480:CACACACCTG:C         CACACACCTGGACAA                  C        NA
#>   1:839515:CTAGACACAC:C CTAGACACAC...CACACACACG                  C        NA
#>                              FILTER
#>                         <character>
#>          1:788418:CAG:C           .
#>       1:789568:TATGGA:T           .
#>   1:790933:CGAATGGAAT:C           .
#>          1:818464:CCT:C           .
#>           1:826577:A:AT           .
#>   1:839480:CACACACCTG:C           .
#>   1:839515:CTAGACACAC:C           .
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome


```
