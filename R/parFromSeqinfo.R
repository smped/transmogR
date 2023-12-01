#' @title Get the PAR-Y Regions From a Seqinfo Object
#'
#' @description Define the Pseudo-Autosomal Regions from a Seqinfo Object
#'
#' @references GenomicRanges object
#'
#' @details
#' Using a seqinfo object based on either hg38, hg19 or their variations,
#' create a GRanges object with the Pseudo-Autosomal Regions from the Y
#' chromosome for that build.
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
#' parFromSeqinfo(sq)
#'
#'
#' @importFrom GenomeInfoDb seqnames
#' @importFrom GenomicRanges GRanges
#' @export
parFromSeqinfo <- function(x, ...){

    ## Retain the seqnames but estimate the genome build.
    stopifnot(is(x, "Seqinfo"))
    hasChr <- any(grepl("^chr",seqnames(x)))
    y <- "Y"
    if (hasChr) y <- "chrY"
    sq <- as.data.frame(x)
    len <- sq[y, "seqlengths"]
    par_df <- data.frame(
        build = c("hg19", "hg38"),
        length = c(59373566, 57227415),
        par1 = c("10001-2781479", "10001-2781479"),
        par2 = c( "59034050-59363566", "56887903-57217415")
    )
    par_df <- subset(par_df, length == len)
    if (nrow(par_df) == 0) stop("Invalid length for ", y)
    rng <- paste0(y, ":", unlist(par_df[c("par1", "par2")]))
    gr <- GRanges(rng, seqinfo = x)
    gr$PAR <- paste0("PAR", seq_len(2))
    gr

}
