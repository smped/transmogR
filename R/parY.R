#' @title Get the PAR-Y Regions From a Seqinfo Object
#'
#' @description Define the Pseudo-Autosomal Regions from a Seqinfo Object
#'
#' @return A GenomicRanges object
#'
#' @details
#' Using a seqinfo object based on either hg38, hg19, CHM13.v2 or their
#' variations, create a GRanges object with the Pseudo-Autosomal Regions from
#' the Y chromosome for that build.
#' The length of the Y chromosome on the seqinfo object is used to determine
#' the correct genome build when passing a Seqinfo object.
#' Otherwise
#'
#' An additional mcols column called PAR will indicate PAR1 and PAR2
#'
#' @param x A Seqinfo object or any of named build. If passing
#' a character vector, [match.arg()] will be used to match the build.
#' @param prefix Optional prefix to place before chromosome names. Can only be
#' NULL, "" or "chr"
#' @param ... Not used
#'
#' @examples
#' library(GenomeInfoDb)
#' sq <- Seqinfo(
#'     seqnames = "chrY", seqlengths = 59373566, genome = "hg19_only_chrY"
#' )
#' parY(sq)
#'
#' ## PAR regions for CHM13 are also available
#' sq <- Seqinfo(
#'     seqnames = "chrY", seqlengths = 62460029, genome = "CHM13"
#' )
#' parY(sq)
#'
#' ## Or just call by name
#' parY("GRCh38", prefix = "chr")
#'
#'
#' @export
#' @name parY
#' @rdname parY-methods
setGeneric("parY", function(x, ...) standardGeneric("parY"))
#' @import GenomicRanges
#' @importFrom GenomeInfoDb seqnames
#' @rdname parY-methods
#' @aliases parY-methods
setMethod(
    "parY", signature = signature(x = "Seqinfo"),
    function(x, ...){

        ## Retain the seqnames but estimate the genome build.
        stopifnot(is(x, "Seqinfo"))
        hasChr <- any(grepl("^chrY",seqnames(x)))
        y <- "Y"
        if (hasChr) y <- "chrY"
        sq <- as.data.frame(x)
        len <- sq[y, "seqlengths"]
        par_df <- data.frame(
            build = c("hg19", "hg38", "chm13v2.0"),
            length = c(59373566, 57227415, 62460029),
            par1 = c("10001-2649520", "10001-2781479", "1-2458320"),
            par2 = c("59034050-59363566", "56887903-57217415", "62122809-62460029")
        )
        par_df <- subset(par_df, length == len)
        if (nrow(par_df) == 0) stop("Invalid length for ", y)
        rng <- paste0(y, ":", unlist(par_df[c("par1", "par2")]))
        gr <- GRanges(rng, seqinfo = x)
        gr$PAR <- paste0("PAR", seq_len(2))
        gr

    }
)
#' @import GenomicRanges
#' @importFrom GenomeInfoDb genome<-
#' @rdname parY-methods
#' @aliases parY-methods
setMethod(
    "parY", signature = signature(x = "character"),
    function(x, prefix = NULL, ...) {

        poss_builds <- c(
            "hg19", "hg38", "chm13v2.0", "GRCh37", "GRCh38", "CHM13"
        )
        x <- match.arg(x, poss_builds)
        if (!is.null(prefix)) stopifnot(prefix %in% c("chr", ""))

        ## Now match the options
        if (x %in% c("hg19", "GRCh37"))
            gr <- paste0(prefix, c("Y:10001-2649520", "Y:59034050-59363566"))

        if (x %in% c("hg38", "GRCh38"))
            gr <- paste0(prefix, c("Y:10001-2781479", "Y:56887903-57217415"))

        if (x %in% c("chr13v2.0", "CHM13"))
            gr <- paste0(prefix, c("Y:1-2458320", "Y:62122809-62460029"))

        gr <- GRanges(gr)
        gr$PAR <- paste0("PAR", seq_len(2))
        genome(gr) <- x
        gr

    }
)
