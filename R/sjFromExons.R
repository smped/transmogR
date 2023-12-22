#' Explore the idea of identifying Splice Junctions from a GTF
library(rtracklayer)
library(plyranges)
gtf <- import.gff(
    system.file("extdata/gencode.v44.subset.gtf.gz", package = "transmogR")
)
gtf <- splitAsList(gtf, gtf$type)
## Reminder that
# - The *donor* site at the 5' end of the intron is usually GG-(cut)-GURAGU where R represents an A or G
#     + Thus we keep two bases of the exon and 6 of the intron
# - The *acceptor* site at the 3' end of the intron is usually a Y-rich-region followed by NCAG-(cut)-G
#     + Thus we keep 1 base of the exon + 4 from the intron
##
sjFromExons <- function(
        x, rank_col = c("exon_number", "exon_rank"),
        tx_col = c("transcript_id", "tx_id"), extra_cols = "all",
        donor_len = 8, acceptor_len = 5, chop = TRUE, ...
){
    stopifnot(is(x, "GRanges"))
    mc <- mcols(x)
    rank_col <- intersect(rank_col, colnames(mc))
    stopifnot(length(rank_col) == 1)
    tx_col <- intersect(tx_col, colnames(mc))
    stopifnot(length(tx_col) == 1)
    if (all(c(extra_cols == "all", !is.null(extra_cols)))) {
        extra_cols <- colnames(mc)
    } else {
        extra_cols <- intersect(extra_cols, colnames(mc))
    }
    extra_cols <- setdiff(extra_cols, c(rank_col, tx_col))

    ## Donor sites can't be the final exon
    df <- data.frame(i = seq_along(x), tx = mc[[tx_col]], rank = mc[[rank_col]])
    df <- dplyr::filter(df, rank < max(rank), .by = tx)
    donor <- x[df$i]
    # Shift downstream
    neg <- strand(donor) == "-"
    donor[neg] <- resize(
        shift(donor[neg], -1 * (donor_len - 2)), donor_len, fix = "end"
    )
    donor[!neg] <- resize(
        shift(donor[!neg], donor_len - 2), donor_len, fix = "end"
    )
    donor$site <- "donor"

    ## Acceptor sites can't be the 1st exon
    acceptor <- x[mc[[rank_col]] != 1]
    ## Shift upstream
    neg <- strand(acceptor) == "-"
    acceptor[neg] <- resize(
        shift(acceptor[neg], acceptor_len - 1), acceptor_len, fix = "start"
    )
    acceptor[!neg] <- resize(
        shift(acceptor[!neg], -1 * (acceptor_len - 1)),
        acceptor_len, fix = "start"
    )
    acceptor$site <- "acceptor"

    sj <- sort(c(donor, acceptor))
    return_cols <- unique(c("site", tx_col, extra_cols, rank_col))
    mcols(sj) <- mcols(sj)[return_cols]
    if (chop) sj <- extraChIPs::chopMC(sj)
    sj

}
sjFromExon(gtf$exon, extra_cols = c(), chop = FALSE)

## Alternatively, just shift everything, then remove
GRangesList(
    donor = gtf$exon %>%
        group_by(transcript_id) %>%
        plyranges::filter(exon_number != max(exon_number)) %>%
        ungroup() %>%
        plyranges::select(gene_id, gene_name, transcript_id, exon_id, exon_number) %>%
        shift_downstream(6) %>%
        resize(8, fix = "end"),
    acceptor = gtf$exon %>%
        subset(exon_number != 1) %>%
        plyranges::select(gene_id, gene_name, transcript_id, exon_id, exon_number) %>%
        shift_upstream(4) %>%
        resize(5, fix = "start")
) %>%
    unlist() %>%
    names_to_column("site") %>%
    plyranges::select(site, gene_id, gene_name, transcript_id, exon_number) %>%
    GenomicRanges::sort() %>%
    extraChIPs::chopMC()

## Checking an ensdb object
library(AnnotationHub)
ah <- AnnotationHub()
ensdb <- ah[["AH113665"]]
ex <- exonsBy(ensdb) %>%
    unlist() %>%
    names_to_column("tx_id") %>%
    subset(seqnames %in% c(1:22, "X", "Y")) %>%
    keepStandardChromosomes() %>%
    sortSeqlevels() %>%
    sort()
ex
# GRanges object with 1813223 ranges and 3 metadata columns:
#     seqnames      ranges strand |           tx_id         exon_id exon_rank
# <Rle>   <IRanges>  <Rle> |     <character>     <character> <integer>
#     [1]        1 11869-12227      + | ENST00000456328 ENSE00002234944         1
# [2]        1 12010-12057      + | ENST00000450305 ENSE00001948541         1
# [3]        1 12179-12227      + | ENST00000450305 ENSE00001671638         2
# [4]        1 12613-12697      + | ENST00000450305 ENSE00001758273         3
# [5]        1 12613-12721      + | ENST00000456328 ENSE00003582793         2
# ...      ...         ...    ... .             ...             ...       ...
# [1813219]  LRG_793 28813-28898      + |       LRG_793t1     LRG_793t1e5         5
# [1813220]  LRG_793 30985-31063      + |       LRG_793t1     LRG_793t1e6         6
# [1813221]  LRG_793 34342-36438      + |       LRG_793t1     LRG_793t1e7         7
# [1813222]   LRG_93   4981-5476      + |        LRG_93t1      LRG_93t1e1         1
# [1813223]   LRG_93 19467-20466      + |        LRG_93t1      LRG_93t1e2         2
# -------
#     seqinfo: 517 sequences (1 circular) from GRCh38 genome

GRangesList(
    donor = ex %>%
        group_by(tx_id) %>%
        plyranges::filter(exon_rank != max(exon_rank)) %>%
        ungroup() %>%
        plyranges::select(tx_id, exon_id, exon_rank) %>%
        shift_downstream(6) %>%
        resize(8, fix = "end"),
    acceptor = ex %>%
        subset(exon_rank != 1) %>%
        plyranges::select(tx_id, exon_id, exon_rank) %>%
        shift_upstream(4) %>%
        resize(5, fix = "start")
) %>%
    unlist() %>%
    names_to_column("site") %>%
    plyranges::select(site, tx_id, exon_rank) %>%
    GenomicRanges::sort() %>%
    extraChIPs::chopMC()

## Looks promising
