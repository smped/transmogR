#' @title Calculate new exon co-ordinates
#'
#' @description Calcluate new exon co-ordinates after including InDels
#'
#' @details
#' Given a set of variants, this will return a set of genomic ranges with
#' updated co-ordinates able to be applied on a variant modified reference
#' genome
#'
#' @param x A GenomicRanges object with co-ordinates needing to be recalculated
#' @param var A set of variants to be incorporated into a reference genome
#' @param alt_col The name of the column with the alternate sequence for each
#' variant
#' @param mc.cores Passed internally to [parallel::mclapply()]
#' @param ... Not used
#'
#' @return
#' GRanges object with co-ordinates shifted according the to the provided
#' variants. The new co-ordinates will be compatible with a variant-modified
#' genome as produced by [genomogrify()] and can be used to extract the
#' sequences associated with the ranges in the modified reference.
#'
#' @examples
#' # Define a 3nt insertion
#' var <- GRanges("seq1:5:*", seqlengths = c(seq1=10), REF = "A", ALT = "AGT")
#' var
#' # A simple GRanges to shift co-ordinates for
#' gr <- GRanges("seq1:1-10:+", seqlengths = c(seq1=10), feature = "feature1")
#' gr
#' # Create shifted co-ordinates based on the provided variants
#' new_gr <- shiftByVar(gr, var)
#' new_gr
#' ## The seqlengths will have been adjusted to account for all variants
#' seqinfo(new_gr)
#'
#' @importFrom GenomeInfoDb seqinfo seqlengths seqlengths<-
#' @importFrom S4Vectors splitAsList
#' @importFrom IRanges IRanges
#' @importFrom parallel mclapply
#' @export
shiftByVar <- function(x, var, alt_col = "ALT", mc.cores = 1, ...) {

    ## The seqinfo object will also require updating after adjusting lengths
    stopifnot(is(x, "GRanges"))
    sq <- seqinfo(x)
    sl <- seqlengths(sq)

    var <- cleanVariants(var, alt_col = alt_col, ...)
    var <- subset(var, seqnames %in% seqlevels(sq))
    var$var_type <- varTypes(var, alt_col = alt_col)
    ## Now keep only the InDels on the required chromosomes
    var <- subset(var, var$var_type != "SNV")
    var_grl <- splitAsList(var, seqnames(var))
    var_grl <- var_grl[vapply(var_grl, length, integer(1)) > 0]
    var_grl <- var_grl[names(var_grl) %in% seqlevels(sq)]

    ## Return the shifts as a list of Rle objects. Doesn't need to be an RleList
    shifts <- lapply(
        var_grl, .createOffsetMap, alt_col = alt_col, seq_lengths = sl
    )
    has_shift <- names(shifts)
    invisible(gc())

    ## Now apply the shifts across each sequence with a variant
    seqlengths(sq) <- NA ## Remove to avoid out-of-bounds errors when shifting
    new_ranges <- GRangesList(
        mclapply(
            has_shift,
            \(i) {
                gr <- subset(x, seqnames == i)
                offset <- shifts[[i]]
                new_start <- start(gr) + as.integer(offset[start(gr)])
                new_ends <- end(gr) + as.integer(offset[end(gr)])
                if (!all(new_ends >= new_start)) {
                    msg <- paste("Unresolvable errors appear present on", i)
                    stop(msg)
                }
                iranges <- IRanges(new_start, new_ends)
                strand <- strand(gr)
                new_gr <- GRanges(i, iranges, strand, seqinfo = sq)
                mcols(new_gr) <- mcols(gr)
                new_gr
            }, mc.cores = mc.cores
        )
    )
    names(new_ranges) <- has_shift
    ## Pull out any ranges from unaffected chromosomes
    old_ranges <- subset(x, !seqnames %in% has_shift)
    old_ranges <- splitAsList(old_ranges, seqnames(old_ranges))
    old_ranges <- old_ranges[!names(old_ranges) %in% has_shift]

    ## Now update the seqlengths
    chr_offsets <- vapply(shifts, \(x) as.integer(x[length(x)]), integer(1))
    sl[has_shift] <- sl[has_shift] + chr_offsets

    ## And return the final object
    out <- unlist(c(old_ranges, new_ranges)[seqlevels(sq)])
    seqlengths(out) <- sl
    out

}

#' @keywords internal
#' @importFrom S4Vectors Rle
#' @importFrom IRanges findOverlaps
#' @importFrom methods slot
.createOffsetMap <- function(gr, alt_col, seq_lengths) {

    ## Handle the deletions using a GPos so individual nucleotides
    ## can be deleted correctly
    dels <- subset(gr, gr$var_type == "Deletion")
    if (length(dels)) {
        dels <- GPos(dels)
        hits <- findOverlaps(dels, gr)
        stopifnot(length(hits) == length(dels)) # Shouldn't happen...
        dels$id <- slot(hits, "to")
        dels$change <- -1 * c(0, diff(start(dels)))
        ## This handles directly neighbouring deletions & also removes the
        ## first position of any deletion
        dels <- as(dels[c(FALSE, diff(dels$id) == 0)], "GRanges")
    }

    ## Insertions are much simpler as a single offset is added for
    ## the individual nt where the insertion occurs
    ins <- subset(gr, gr$var_type == "Insertion")
    ins$change <- nchar(mcols(ins)[[alt_col]]) - nchar(ins$REF)

    ## Now create the map
    map <- sort(c(ins, dels))
    seq_name <- unique(seqnames(gr))
    n <- seq_lengths[[seq_name]]
    if (is.na(n) | n == 0)
        stop("The seqinfo object needs to be defined with seqlengths")
    change <- Rle(0, n)
    change[start(map)] <- map$change
    ## Return an Rle of change as the cumulative shift
    cumsum(change)
}


