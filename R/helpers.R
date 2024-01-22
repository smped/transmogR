#' @keywords internal
#' @importFrom Biostrings IUPAC_CODE_MAP
#' @importFrom S4Vectors mcols mcols<-
.checkAlts <- function(var, alt_col) {

    alt_col <- match.arg(alt_col, colnames(mcols(var)))

    ## Deal with NAs
    alts <- mcols(var)[[alt_col]]
    if (is(alts, "XStringSetList")) alts <- as.character(unlist(alts))
    is_na <- is.na(alts)
    if (any(is_na)) message("NA values found in ", alt_col, ". Removing...")
    var <- var[!is_na]
    if (all(is_na) & length(is_na) > 0) stop("All NA values in ", alt_col)

    ## Check Non-IUPAC alt alleles
    iupac <- paste(names(IUPAC_CODE_MAP), collapse = "")
    pat <- paste0("^[", iupac, "]+$")
    alt_error <- !grepl(pat, alts)
    if (any(alt_error)) stop("Non-IUPAC alleles detected")
    mcols(var)[[alt_col]] <- alts[!is_na]
    var
}

#' @keywords internal
#' @import GenomicRanges
#' @importClassesFrom VariantAnnotation ScanVcfParam
#' @importFrom VariantAnnotation readVcf ScanVcfParam vcfWhich<-
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom SummarizedExperiment rowRanges
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

#' @keywords internal
#' @importFrom ComplexUpset intersection_size
.makeIntersectionArgs <- function(x){
    stopifnot(is(x, "list"))
    args <- as.list(formals(intersection_size))
    args <- args[names(args) != "..."]
    cmn <- intersect(names(x), names(args))
    novel <- setdiff(names(x), names(args))
    if (length(cmn) > 0) args[cmn] <- x[cmn]
    if (length(novel) > 0) args <- c(args, x[novel])
    args
}

#' @importFrom S4Vectors mcols mcols<-
#' @importFrom InteractionSet GInteractions
#' @keywords internal
.giFromSj <- function(sj, tx_col, rank_col) {

    stopifnot("site" %in% colnames(mcols(sj)))
    stopifnot(all(sj$site %in% c("donor", "acceptor")))
    stopifnot(is.numeric(mcols(sj)[[rank_col]]))

    site <- c() ## R CMD check error avoidance
    dnr <- subset(sj, site == "donor")
    dnr$sj <- mcols(dnr)[[rank_col]]
    dnr_ids <- paste(mcols(dnr)[[tx_col]], dnr$sj)

    acc <- subset(sj, site == "acceptor")
    acc$sj <- mcols(acc)[[rank_col]] - 1
    acc_ids <- paste(mcols(acc)[[tx_col]], acc$sj)

    dnr_to_acc <- match(dnr_ids, acc_ids)
    acc_to_dnr <- match(acc_ids, dnr_ids)
    if (any(is.na(dnr_to_acc))) stop("NA values when mapping junctions")

    cols <- setdiff(colnames(mcols(dnr)), c("site", rank_col))
    gi <- GInteractions(
        anchor1 = granges(dnr)[dnr_to_acc],
        anchor2 = granges(acc)[acc_to_dnr]
    )
    ## Single columns tend to misbehave a little when adding mcols
    if (length(cols) > 1) {
        mcols(gi) <- mcols(dnr)[dnr_to_acc, cols]
    } else {
        mcols(gi)[[cols]] <- mcols(dnr)[dnr_to_acc, cols]
    }
    gi
}
