#' transmogR: Create a variant-modified reference transcriptome
#'
#' The package `transmogR` has been designed for creation of a variant-modified
#' reference transcriptome
#'
#' The package `transmogR` provides two primary functions for modifying complete
#' transcriptomes or genomes:
#'
#' * [transmogrify()] for incorporating the supplied variants into
#' transcriptomic sequences, and
#' * [genomogrify()] for incorporating the supplied variants into genomic
#' sequences, ideally to be passed as decoy sequences to a tool such as `salmon`.
#'
#' The main functions rely on lower-level functions such as:
#'
#' * [owl()] which over-writes letters (i.e. SNPs) within a sequence, and
#' * [indelcator()] which incorporates InDels into an individual sequence
#'
#' Additional utility functions are provided which allow characterisation and
#' exploration of any set of variants:
#'
#' * [overlapsByVar()] counts the variants which overlap sets of GenomicRanges,
#' first splitting the variants into SNV, Insertions and Deletions
#' * [parY()] returns the pseudo-autosomal regions for a chosen genome build as
#' a GenomicRanges object
#' * [upsetVarByCol()] produces an UpSet plot counting how many unique IDs are
#' impacted by a set o variants. IDs can represent any column in the supplied
#' ranges, such as gene_id or transcript_id
#' * [varTypes()] classifies a set of variants into SNV, Insertions of Deletions
#'
#' @author
#' Stevie Pederson
#'
#' @keywords internal
"_PACKAGE"
