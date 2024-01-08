#' @title Obtain Splice-Junctions from Exons and Transcripts
#'
#' @description
#' Using GRanges defining exons and transcripts, find the splice-junctions
#'
#' @details
#' A canonical splice junction consists of a donor site and an acceptor site at
#' each end of an intron, with a branching site somewhere wthin the intron.
#' Canonical donor sites are 8nt long with the the first two bases being exonic
#' and the next 6 being derived form intronic sequences.
#' Canonical acceptor sites are 5nt long with the first four bases being
#' intronic and the final base being the first base of the next exon.
#'
#' This functions uses each set of exons within a transcript to identify both
#' donor and acceptor sites. Branch sites are not identified.
#'
#' @return
#' A GRanges object with requested columns, and an additional column, 'site',
#' annotating each region as a donor or acceptor site.
#'
#' Alternatively, by specifying as = "GInteractions", the junctions can be
#' returned with each splice junction annotated as a GenomicInteraction.
#' This can make the set of junctions easier to interpret for a given transcript.
#'
#' @param x GRanges object with exons and transcripts. A column indicating the
#' position (or rank) of each exon within the transcript must be included.
#' @param rank_col The column containing the position of each exons within the
#' transcript
#' @param tx_col The column containing unique transcript-level identifiers
#' @param extra_cols Can be a vector of column names to return beyond rank_col
#' and tx_col. By default all columns are returned (extra_cols = "all").
#' @param don_len,acc_len Length of donor and acceptor sites respectively
#' @param as Return as a set of GenomicRanges, or with each splice junction
#' annotated as a GenomicInteraction
#' @param ... Not used
#'
#' @examples
#' library(rtracklayer)
#' gtf_cols <- c(
#'   "transcript_id", "transcript_name", "gene_id", "gene_name", "exon_number"
#' )
#' gtf <- import.gff(
#'    system.file("extdata/gencode.v44.subset.gtf.gz", package = "transmogR"),
#'    feature.type = "exon", colnames = gtf_cols
#' )
#' sj <- sjFromExons(gtf)
#' sj
#'
#' ## Or to simplify shared splice junctions across multiple transcripts
#' library(extraChIPs, quietly = TRUE)
#' chopMC(sj)
#'
#' ## Splice Junctions can also be returned as a GInteractions object with
#' ## anchorOne as the donor & anchorTwo as the acceptor sites
#' sjFromExons(gtf, as = "GInteractions")
#'
#'
#' @importFrom S4Vectors mcols 'mcols<-'
#' @importFrom methods is
#' @importFrom stats aggregate
#' @export
sjFromExons <- function(
        x, rank_col = c("exon_number", "exon_rank"),
        tx_col = c("transcript_id", "tx_id"), extra_cols = "all",
        don_len = 8, acc_len = 5, as = c("GRanges", "GInteractions"), ...
){
    stopifnot(is(x, "GRanges"))
    as <- match.arg(as)
    mc <- mcols(x)
    rank_col <- intersect(rank_col, colnames(mc))
    stopifnot(length(rank_col) == 1)
    tx_col <- intersect(tx_col, colnames(mc))
    stopifnot(length(tx_col) == 1)
    if (any(strand(x) == "*")) stop("Exons and transcripts must be stranded")
    if (all(c(extra_cols == "all", !is.null(extra_cols)))) {
        extra_cols <- colnames(mc)
    } else {
        extra_cols <- intersect(extra_cols, colnames(mc))
    }
    extra_cols <- setdiff(extra_cols, c(rank_col, tx_col))

    ## Donor sites can't be the final exon
    df <- data.frame(i = seq_along(x), tx = mc[[tx_col]], ex = mc[[rank_col]])
    df$ex <- as.integer(df$ex)
    fm <- as.formula("ex ~ tx")
    max_ex <- aggregate(fm, data = df, FUN = max)
    omit_tx_ex <- paste(max_ex$tx, max_ex$ex, sep = "_")
    all_tx_ex <- paste(df$tx, df$ex, sep = "_")
    df <- df[!all_tx_ex %in% omit_tx_ex,]
    dnr <- x[df$i]
    # Shift downstream
    neg <- strand(dnr) == "-"
    dnr[neg] <- resize(
        shift(dnr[neg], -1 * (don_len - 2)), don_len, fix = "end"
    )
    dnr[!neg] <- resize(shift(dnr[!neg], don_len - 2), don_len, fix = "end")
    dnr$site <- "donor"

    ## Acceptor sites can't be the 1st exon
    acc <- x[mc[[rank_col]] != 1]
    ## Shift upstream
    neg <- strand(acc) == "-"
    acc[neg] <- resize(shift(acc[neg], acc_len - 1), acc_len, fix = "start")
    acc[!neg] <- resize(
        shift(acc[!neg], -1 * (acc_len - 1)), acc_len, fix = "start"
    )
    acc$site <- "acceptor"

    sj <- sort(c(dnr, acc))
    return_cols <- unique(c("site", tx_col, extra_cols, rank_col))
    mcols(sj)[[rank_col]] <- as.integer(mcols(sj)[[rank_col]])
    mcols(sj) <- mcols(sj)[return_cols]
    if (as == "GInteractions") sj <- .giFromSj(sj, tx_col, rank_col)
    sj

}

