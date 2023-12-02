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
