#' @title Count overlaps by variant type
#'
#' @description
#' Count how many variants overlap each element of a GRangesList
#'
#' @details
#' Taking any GRangesList, count how many of each variant type overlap a region
#' within each element.
#'
#' @return A data.frame
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
#' @importFrom S4Vectors mcols splitAsList
#' @importFrom IRanges overlapsAny
#' @export
overlapsByVar <- function(x, var, alt_col = "ALT", ...) {
    stopifnot(is(x, "GRangesList"))
    stopifnot(is(var, "GRanges"))
    alt_col <- match.arg(alt_col, colnames(mcols(var)))
    type <- calvInDel(var, alt_col)
    var_list <- splitAsList(var, type)
    df_list <- lapply(
        x,
        \(i) {
            res <- lapply(var_list, overlapsAny, i)
            as.data.frame(lapply(res, sum))
        }
    )
    df <- do.call("rbind", df_list)
    df$Total <- rowSums(df, ...)
    df
}
