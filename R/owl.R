#' @title OverWrite Letters in an XStringSet
#'
#' @description
#' OverWrite Letters (e.g. SNPs) in an XStringSet
#'
#' @param seq A BSgenome, DNAStringSet, RNAStringSet or other XStringSet.
#' @param snps A GRanges object with SNP positions and a column containing
#' the alternate allele
#' @param alt_col Column name in the mcols element of `snps` containing the
#' alternate allele
#' @param ... Passed to [Biostrings::replaceLetterAt()]
#'
#' @details
#' If providing a BSgenome object, this will first be coerced to a DNAStringSet
#' which can be time consuming
#'
#' @return An object of the same class as the original object, but with SNPs
#' inserted at the supplied positions
#'
#' @examples
#' library(GenomicRanges)
#' seq <- DNAStringSet(c(chr1 = "AAGC"))
#' snps <- GRanges("chr1:2")
#' snps$ALT <- "G"
#' snps
#' seq
#' owl(seq, snps)
#'
#' @import Biostrings
#' @importFrom S4Vectors mcols
#' @importFrom methods is
#' @importFrom GenomeInfoDb seqinfo 'seqinfo<-' seqlevels seqnames 'seqlevels<-'
#' @importFrom BSgenome getSeq
#'
#' @export
#' @name owl
#' @rdname owl-methods
setGeneric(
  "owl", function(seq, snps, ...){standardGeneric("owl")}
    ## A more boring name might be subSNP...
)
#' @import Biostrings
#' @importFrom S4Vectors mcols
#' @importFrom methods is
#' @importFrom GenomeInfoDb seqinfo 'seqinfo<-' seqlevels seqnames 'seqlevels<-'
#' @rdname owl-methods
#' @aliases owl
#' @export
setMethod(
  "owl",
  signature = signature(seq = "XStringSet", snps = "GRanges"),
  function(seq, snps, alt_col = "ALT", ...) {

    ## Check the SNPs
    stopifnot(is(snps, "GenomicRanges"))
    stopifnot(all(width(snps) == 1)) # Must be single positions
    alt_col <- match.arg(alt_col, colnames(mcols(snps)))
    stopifnot(all(nchar(mcols(snps)[[alt_col]]) == 1))

    ## Check compatible seqinfo
    seq_sq <- seqinfo(seq)
    snp_sq <- seqinfo(snps)
    stopifnot(all(seqlevels(snp_sq) %in% seqlevels(seq_sq))) # seqlevels
    seqlevels(snps) <- seqlevels(seq_sq)
    seqinfo(snps) <- seq_sq

    ## Find sequences with SNPs
    seqs_wth_snps <- as.character(unique(seqnames(snps)))
    new_seq <- lapply(
      seqs_wth_snps,
      function(x) {
        temp <- subset(snps, seqnames == x)
        replaceLetterAt(seq[[x]], at = start(temp), mcols(temp)[[alt_col]], ...)
      }
    )
    names(new_seq) <- seqs_wth_snps
    seq_len <- vapply(seq, length, integer(1))
    seq[seqs_wth_snps] <- new_seq
    ## Check lengths aren't messed up then return
    stopifnot(all(seq_len == vapply(seq, length, integer(1))))
    seq

  }
)
#' @importFrom BSgenome getSeq
#' @rdname owl-methods
#' @aliases owl
#' @export
setMethod(
  "owl",
  signature = signature(seq = "BSgenome", snps = "GRanges"),
  function(seq, snps, alt_col = "ALT", ...) {

    ## Setup the sequence info
    message("Extracting sequences as a DNAStringSet...", appendLF = FALSE)
    seq <- getSeq(seq)
    message("done")
    owl(seq, snps, alt_col = alt_col, ...)
  }
)
