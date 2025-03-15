#' @title Check provided variants for compatibility
#'
#' @description Check and cleanup a set of variants for downstream compatibility
#'
#' @details
#' This function checks a set of variants for the expected structure which is
#' required by all downstream functions in transmogR.
#' The primary change to the data structure is that both REF and ALT columns
#' will be set as simple character vectors.
#'
#' Given the complicated variant calls that can often be produced by variant
#' callers, additional checks performed will be to ensure that:
#'
#' - there are no overlapping variants
#' - SNPs are all single nucleotides and not longer strings
#' - SNPs are bi-allelic
#' - Insertions contain a single nucleotide in the REF column
#' - Deletions contain a single nucleotide in the ALT column
#' - No missing values are present
#' - All ALT/REF nucleotides conform to IUPAC codings
#'
#' @param var GRanges object containing the variants, or a
#' [VariantAnnotation::VcfFile]
#' @param ref_col,alt_col Column names corresponding to the reference and
#' alternate alleles
#' @param ol_vars Error handling for any overlapping variants. Can take values in
#' c("fail", "none", "first", "last", "longest", "shortest").
#' Default is set to fail, with additional options to drop all overlapping
#' variants ('none'), select by genomic position ('first', 'last'), or select
#' by the scale of change to the genome ('longest', 'shortest')
#' @param which Passed to [VariantAnnotation::ScanVcfParam] if working with a
#' VcfFile
#' @param ... Not used
#'
#' @return
#' GRanges object with any incompatible variants removed, or an error produced.
#' The mcols will contain the columns REF and ALT, unless otherwise specified,
#' as character vectors
#'
#' @examples
#' # Any conflicting variants will be removed
#' var <- GRanges(c("chr10:114468420-114468422", "chr10:114468422"))
#' var$REF <- c("GCC", "C")
#' var$ALT <- c("G", "CTAT")
#' var
#' ## Taken from the 1000GP, the first variant would delete the C at 114468422
#' ## whilst the second variant begins an insertion at this position.
#' ## These are clearly conflicting. The default value for ol_vars is to fail
#' ## with an error (ol_vars = "fail"). However, both can be removed by setting
#' ## ol_vars = "none". A warning will always be produced.
#' cleanVariants(var, ol_vars = "none")
#' ## Or the longest can be retained, along with multiple other options
#' cleanVariants(var, ol_vars = "longest")
#'
#' @export
#' @name cleanVariants
#' @rdname cleanVariants-methods
setGeneric("cleanVariants", function(var, ...) standardGeneric("cleanVariants"))
#' @export
#' @rdname cleanVariants-methods
#' @aliases cleanVariants-methods
setMethod(
    "cleanVariants",
    signature = signature(var = "GRanges"),
    function(
        var, ol_vars = "fail", ref_col = "REF", alt_col = "ALT", ...
    ){
        var <- .checkAlts(var, alt_col, ref_col, ol_vars)
        .checkOverlapVars(var, ol_vars)
    }
)
#' @importClassesFrom VariantAnnotation VcfFile
#' @export
#' @rdname cleanVariants-methods
#' @aliases cleanVariants-methods
setMethod(
    "cleanVariants",
    signature = signature(var = "VcfFile"),
    function(
        var, ol_vars = "fail", which, ref_col = "REF", alt_col = "ALT", ...
    ) {
        var <- .parseVariants(var, alt_col, which) ## Check args here
        cleanVariants(var, ol_vars, ref_col, alt_col, ...)
    }
)
