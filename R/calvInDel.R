#' @title Identify SNVs, Insertions and Deletions
#'
#' @description
#' Identify SNVs, Insertions and Deletions within a GRanges object
#'
#' @details
#' Using the width of the reference and alternate alleles, classify each
#' range as an SNV, Insertion or Deletion.
#'
#' - SNVs are expected to have REF & ALT widths of 1
#' - Insertions are expected to have ALT longer than REF
#' - Deletions are expected to have ALT shorter than REF
#'
#' These are relatively permissive criteria
#'
#' @return Character vector
#'
#' @param x GenomicRanges object
#' @param alt_col Name of the column with mcols(x) which contains the alternate
#' allele. Can be an XStringSetList, XStringSet or character
#' @param ... Not used
#'
#' @examples
#' # Load the example VCF and classify ranges
#' library(VariantAnnotation)
#' f <- system.file("extdata/1000GP_subset.vcf.gz", package = "transmogR")
#' vcf <- readVcf(f)
#' gr <- rowRanges(vcf)
#' type <- calvInDel(gr)
#' table(type)
#' gr[type != "SNV"]
#'
#'
#'
#' @importFrom S4Vectors mcols
#' @importFrom methods is
#' @importFrom IRanges width
#' @export
calvInDel <- function(x, alt_col = "ALT", ...){
    ## x should be a GRanges object with variants
    stopifnot(is(x, "GRanges"))
    stopifnot(alt_col %in% colnames(mcols(x)))

    ## Basic info
    w <- width(x)
    n <- length(x)

    ## Check the ALTs
    alts <- mcols(x)[[alt_col]]
    if (is(alts, "XStringSetList")) alts <- unlist(alts)
    stopifnot(is(alts, "XStringSet") | is(alts, "character"))
    alt_width <- nchar(alts)

    ## Classify each range
    snv <- w == 1 & alt_width == 1
    ins <- alt_width > w
    dels <- w > alt_width
    ## All should add up to the total
    if (sum(snv + ins + dels) != n)
        stop("Some variants were unable to be uniquely identified")

    ## Return output
    out <- character(n)
    out[snv] <- "SNV"
    out[ins] <- "Insertion"
    out[dels] <- "Deletion"
    out
}
