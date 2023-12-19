#' @keywords internal
#' @importFrom Biostrings IUPAC_CODE_MAP
#' @importFrom S4Vectors mcols
.checkAlts <- function(var, alt_col) {

    alt_col <- match.arg(alt_col, colnames(mcols(var)))

    ## Deal with NAs
    is_na <- is.na(mcols(var)[[alt_col]])
    if (any(is_na)) message("NA values found in ", alt_col, ". Removing...")
    var <- var[!is_na]

    ## Check Non-IUPAC alt alleles
    iupac <- paste(names(IUPAC_CODE_MAP), collapse = "")
    pat <- paste0("^[", iupac, "]+$")
    alt_error <- !grepl(pat, mcols(var)[[alt_col]])
    if (any(alt_error)) stop("Non-IUPAC alleles detected")
    var
}

#' @keywords internal
#' @importClassesFrom VariantAnnotation ScanVcfParam
#' @importFrom VariantAnnotation readVcf ScanVcfParam 'vcfWhich<-'
#' @importFrom S4Vectors mcols 'mcols<-'
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb seqinfo
#' @importFrom methods is
.parseVariants <- function(f, alt_col, which, ...){
    param <- ScanVcfParam(fixed = alt_col, info = NA, ...)
    if (!missing(which)) {
        stopifnot(is(which, "GRanges"))
        vcfWhich(param) <- which
    }
    vcf <-  readVcf(f, param = param)
    gr <- rowRanges(vcf)
    mc_names <- c("REF", alt_col)
    mcols(gr) <- mcols(gr)[mc_names]
    mcols(gr) <- lapply(
        mcols(gr), \(x) {
            if (is(x, "XStringSetList")) x <- unlist(x)
            as.character(x)
        }
    )
    gr
}
