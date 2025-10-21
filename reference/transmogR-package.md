# transmogR: Create a variant-modified reference transcriptome

The package `transmogR` has been designed for creation of a
variant-modified reference transcriptome

## Details

The package `transmogR` provides two primary functions for modifying
complete transcriptomes or genomes:

- [`transmogrify()`](https://smped.github.io/transmogR/reference/transmogrify-methods.md)
  for incorporating the supplied variants into transcriptomic sequences,
  and

- [`genomogrify()`](https://smped.github.io/transmogR/reference/genomogrify-methods.md)
  for incorporating the supplied variants into genomic sequences,
  ideally to be passed as decoy sequences to a tool such as `salmon`.

The main functions rely on lower-level functions such as:

- [`owl()`](https://smped.github.io/transmogR/reference/owl-methods.md)
  which over-writes letters (i.e. SNPs) within a sequence, and

- [`indelcator()`](https://smped.github.io/transmogR/reference/indelcator-methods.md)
  which incorporates InDels into an individual sequence

Additional utility functions are provided which allow characterisation
and exploration of any set of variants:

- [`overlapsByVar()`](https://smped.github.io/transmogR/reference/overlapsByVar-methods.md)
  counts the variants which overlap sets of GenomicRanges, first
  splitting the variants into SNV, Insertions and Deletions

- [`parY()`](https://smped.github.io/transmogR/reference/parY-methods.md)
  returns the pseudo-autosomal regions for a chosen genome build as a
  GenomicRanges object

- [`upsetVarByCol()`](https://smped.github.io/transmogR/reference/upsetVarByCol.md)
  produces an UpSet plot counting how many unique IDs are impacted by a
  set o variants. IDs can represent any column in the supplied ranges,
  such as gene_id or transcript_id

- [`varTypes()`](https://smped.github.io/transmogR/reference/varTypes.md)
  classifies a set of variants into SNV, Insertions of Deletions

## See also

Useful links:

- <https://github.com/smped/transmogR>

- Report bugs at <https://github.com/smped/transmogR/issues>

## Author

Stevie Pederson
