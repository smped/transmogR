# Creation of example gtf from Gencode v44

The original file was obtained from https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz
The subset was then created using

```r
library(VariantAnnotation)
library(rtracklayer)
library(magrittr)
vcf <- VcfFile(
    system.file("extdata/1000GP_subset.vcf.gz", package = "transmogR")
)
gr <- transmogR:::.parseVariants(vcf, "ALT")
gtf <- import.gff("gencode.v44.primary_assembly.annotation.gtf.gz")
ids <- unique(subsetByOverlaps(gtf, gr)$gene_id)
subset(gtf, gene_id %in% ids) %>%
  export.gff("inst/extdata/gencode.v44.subset.gtf", format = "gtf")
```

Then moving to `bash`

```bash
gzip inst/extdata/gencode.v44.subset.gtf
```
