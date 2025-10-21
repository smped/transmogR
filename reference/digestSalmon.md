# Parse the output from salmon

Parse transcript counts and additional data from salmon

Calculate the overdispersions from a set of paths without parsing any
counts

## Usage

``` r
digestSalmon(
  paths,
  max_sets = 2L,
  aux_dir = "aux_info",
  name_fun = basename,
  verbose = TRUE,
  extra_assays = NULL,
  max_boot = Inf,
  ...
)

overdispFromBoots(paths, n_boot, .ids)
```

## Arguments

- paths:

  Vector of file paths to directories containing salmon results

- max_sets:

  The maximum number of indexes permitted

- aux_dir:

  Subdirectory where bootstraps and meta_info.json are stored

- name_fun:

  Function applied to paths to provide colnames in the returned object.
  Set to NULL or c() to disable.

- verbose:

  Print progress messages

- extra_assays:

  Can take values in c("TPM", "effectiveLength", "length") to optionally
  request TPM, effectiveLength or length as assays. Including the length
  assay is intended for the use case of personalised transcriptomes
  where transcript lengths may no longer be uniform across samples. None
  will be returned by default

- max_boot:

  The maximum number of bootstraps to use. Setting this to zero will
  ignore all bootstraps and the scaledCounts assay will not be included
  in the returned object

- ...:

  Not used

- n_boot:

  The number of bootstraps

- .ids:

  Vector of transcript IDs which match the bootstrap values. Will be
  parsed from paths if not provided, although this adds time

## Value

A SummarizedExperiment object containing assays for counts and
scaledCounts. The scaledCounts assay contains counts divided by
overdispersions. rowData in the returned object will also include
transcript-lengths along with the overdispersion estimates used to
return the scaled counts. TPM, effectiveLength and length can be
returned as additional assays by specifying one or more of these in the
extra_assays argument

`overdispFromBoots` returns a numeric vector

## Details

This function is based heavily on
[`edgeR::catchSalmon()`](https://rdrr.io/pkg/edgeR/man/catchSalmon.html)
however, there are some important differences:

1.  A SummarizedExperiment object is returned

2.  Differing numbers of transcripts are allowed between samples

The second point is intended for the scenario where some samples may
have been aligned to a full reference, with remaining samples aligned to
a partially masked reference (e.g. chrY). This will lead to differing
numbers of transcripts within each salmon index, however, common
estimates of overdispersions are required for scaling transcript-level
counts. By default, the function will error if \>2 different sets of
transcripts are detected, however this can be modified using the
max_sets argument.

This greater flexibility also requires more stringent checking and, as
such, for smaller datasets, digestSalmon may be slower that the edgeR
function.

The SummarizedExperiment object returned may also contain multiple
assays, as described elsewhere on this page

This follows the methods of Baldoni, et al. (2024). Dividing out
quantification uncertainty allows efficient assessment of differential
transcript expression with edgeR. Nucleic Acids Research, 52(3), e13.
https://doi.org/10.1093/nar/gkad1167

## Examples

``` r
## Provide the path to the parent directories which contains each
## quant.sf file
ex_path <- system.file("extdata/salmon_test", package = "transmogR")
se <- digestSalmon(ex_path, extra_assays = "TPM", verbose = FALSE)
se
#> class: SummarizedExperiment 
#> dim: 2 1 
#> metadata(2): resampleType n_boot
#> assays(3): counts scaledCounts TPM
#> rownames(2): ENST00000000233.10 ENST00000000412.8
#> rowData names(2): overdispersion length
#> colnames(1): salmon_test
#> colData names(2): totals n_trans

ex_path <- system.file("extdata/salmon_test", package = "transmogR")
overdispFromBoots(ex_path, 10)
#> [1]  80.7445 179.9536
```
