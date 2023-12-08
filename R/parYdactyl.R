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
#' the correct genome build
#'
#' An additional mcols column called PAR will indicate PAR1 and PAR2
#'
#' @param x A Seqinfo object
#' @param ... Not used
#'
#' @examples
#' library(GenomeInfoDb)
#' sq <- Seqinfo(
#'   seqnames = "chrY", seqlengths = 59373566, genome = "hg19_only_chrY"
#' )
#' parYdactyl(sq)
#'
#'
#' @importFrom GenomeInfoDb seqnames
#' @importFrom GenomicRanges GRanges
#' @export
parYdactyl <- function(x, ...){

    # Original name was parFromSeqinfo

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
        par1 = c("10001-2781479", "10001-2781479", "0-2458320"),
        par2 = c("59034050-59363566", "56887903-57217415", "62122809-62460029")
    )
    par_df <- subset(par_df, length == len)
    if (nrow(par_df) == 0) stop("Invalid length for ", y)
    rng <- paste0(y, ":", unlist(par_df[c("par1", "par2")]))
    gr <- GRanges(rng, seqinfo = x)
    gr$PAR <- paste0("PAR", seq_len(2))
    gr

}
