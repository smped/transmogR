#' @keywords internal
#' @importFrom methods slot
.checkOverlapVars <- function(
        var,
        ol_vars = c("fail", "none", "first", "last", "longest", "shortest")
){

    ol_vars <- match.arg(ol_vars)
    reduced_gr <- reduce(var, ignore.strand = TRUE, min.gapwidth = 0L)
    ol <- findOverlaps(reduced_gr, var)
    dups <- duplicated(slot(ol, "from"))
    reduced_multi <- unique(slot(ol, "from")[dups])
    self_ol <- subset(ol, queryHits %in% reduced_multi)

    if (length(self_ol) > 0) {
        msg <- paste(
            length(self_ol),
            "pairs of overlapping loci found and cannot be incorporated into a modified reference."
        )
        if (ol_vars == "fail") stop(msg)
        if (ol_vars == "none") {
            msg <- c(msg, "\nAll overlapping variants will be removed")
            keep <- NULL
        }

        ## From this point a selection will be returned so setup a list
        ol_list <- split(slot(self_ol, "to"), slot(self_ol, "from"))
        if (ol_vars == "first") {
            msg <- c(msg, "\nThe first overlapping locus by genomic position will be retained")
            keep <- vapply(ol_list, \(x) x[1], integer(1))
        }
        if (ol_vars == "last") {
            msg <- c(msg, "\nThe last overlapping locus by genomic position will be retained")
            keep <- vapply(ol_list, \(x) x[length(x)], integer(1))
        }

        ## The remaining two options will require a list with lengths
        if (ol_vars == "longest") {
            msg <- c(msg, "\nThe longest overlapping locus by genomic change will be retained")
            keep <- lapply(
                ol_list,
                \(x) x[which.max(abs(nchar(var[x]$REF) - nchar(var[x]$ALT)))]
            )
            keep <- unlist(keep)
        }
        if (ol_vars == "shortest") {
            msg <- c(msg, "\nThe shortest overlapping locus by genomic change will be retained")
            keep <- lapply(
                ol_list,
                \(x) x[which.min(abs(nchar(var[x]$REF) - nchar(var[x]$ALT)))]
            )
            keep <- unlist(keep)
        }

        discard <- unique(setdiff(slot(self_ol, "to"), keep))
        var <- var[-discard]
        warning(msg)
    }
    var
}

#' @keywords internal
#' @importFrom Biostrings IUPAC_CODE_MAP
#' @importFrom S4Vectors mcols mcols<-
.checkAlts <- function(var, alt_col, ref_col = "REF", ol_vars) {

    alt_col <- match.arg(alt_col, colnames(mcols(var)))
    alts <- mcols(var)[[alt_col]]
    if (is(alts, "XStringSetList")) alts <- as.character(unlist(alts))

    ## Check duplicate loci
    err <- c()
    dup <- duplicated(var)
    if (any(dup))
        err <- c(
            err, "Duplicate variant loci found. Please choose which one to use"
        )

    ## Multi allelic sites
    if (length(alts) != length(var)) {
        err <- c(
            err, "Alternate alleles have been specified with multiple values"
        )
    }

    ## Deal with NAs
    is_na <- is.na(alts)
    if (any(is_na)) err <- c(err, "NA values found in alternate alleles")

    ## Empty deletions
    if (any(nchar(alts[!is_na]) == 0)) {
        err <- c(
            err,
            "Please set Deletions so that width(REF) > 1 and width(ALT) > 0"
        )
    }

    ## Non-IUPAC alt alleles
    iupac <- paste(names(IUPAC_CODE_MAP), collapse = "")
    pat <- paste0("^[", iupac, "]+$")
    alt_error <- !grepl(pat, alts[!is_na])
    if (any(alt_error))
        err <- c(err, "Non-IUPAC alleles detected (may include empty ALTs)")
    ## Report errors
    if (length(err) > 0) {
        err <- paste(err, collapse = "\n\t")
        stop("\nFound:\n\t", err)
    }

    mcols(var)[[alt_col]] <- alts[!is_na]
    var

}

#' @keywords internal
#' @import GenomicRanges
#' @importClassesFrom VariantAnnotation ScanVcfParam
#' @importFrom VariantAnnotation readVcf ScanVcfParam vcfWhich<-
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom Seqinfo seqinfo
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
.makeIntersectionArgs <- function(x){
    stopifnot(is(x, "list"))
    args <- as.list(formals(ComplexUpset::intersection_size))
    args <- args[names(args) != "..."]
    cmn <- intersect(names(x), names(args))
    novel <- setdiff(names(x), names(args))
    if (length(cmn) > 0) args[cmn] <- x[cmn]
    if (length(novel) > 0) args <- c(args, x[novel])
    args
}

#' @importFrom S4Vectors mcols mcols<-
#' @keywords internal
.giFromSj <- function(sj, tx_col, rank_col) {

    if (!requireNamespace('InteractionSet', quietly = TRUE))
        stop("Please install 'InteractionSet' to return a GInteractions object.")

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
    gi <- InteractionSet::GInteractions(
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




