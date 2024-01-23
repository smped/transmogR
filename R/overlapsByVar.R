#' @title Count overlaps by variant type
#'
#' @description
#' Count how many variants of each type overlap ranges
#'
#' @details
#' Taking any GRanges or GRangesList, count how many of each variant type
#' overlap a region.
#'
#' @return A vector or matrix
#'
#' @param x A GRangesList with features of interest
#' @param var A Granges object with variants of interest
#' @param alt_col The column within mcols(var) which contains the alternate
#' allele
#' @param ... Passed to \link{rowSums}
#'
#' @examples
#' library(rtracklayer)
#' library(VariantAnnotation)
#' gtf <- import.gff(
#'     system.file("extdata/gencode.v44.subset.gtf.gz", package = "transmogR")
#' )
#' grl <- splitAsList(gtf, gtf$type)
#' vcf <- system.file("extdata/1000GP_subset.vcf.gz", package = "transmogR")
#' var <- rowRanges(readVcf(vcf, param = ScanVcfParam(fixed = "ALT")))
#' overlapsByVar(grl, var)
#'
#' @export
#' @name overlapsByVar
#' @rdname overlapsByVar-methods
setGeneric(
    "overlapsByVar", function(x, var, ...){standardGeneric("overlapsByVar")}
)
#'
#' @importFrom IRanges overlapsAny
#' @export
#' @rdname overlapsByVar-methods
#' @aliases overlapsByVar-methods
setMethod(
    "overlapsByVar", signature = signature(x = "GRangesList", var = "GRanges"),
    function(x, var, alt_col = "ALT", ...) {
        l <- lapply(x, overlapsByVar, var, alt_col = alt_col)
        mat <- do.call("rbind", l)
        cbind(mat, Total = rowSums(mat))
    }
)
#' @importFrom IRanges overlapsAny
#' @export
#' @rdname overlapsByVar-methods
#' @aliases overlapsByVar-methods
setMethod(
    "overlapsByVar", signature = signature(x = "GRanges", var = "GRanges"),
    function(x, var, alt_col = "ALT", ...) {
        ol <- overlapsAny(var, x)
        type <- varTypes(var, alt_col)
        vapply(split(ol, type), sum, integer(1))
    }
)


