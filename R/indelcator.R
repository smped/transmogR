#' @title Substitute InDels into one or more sequences
#'
#' @description
#' Modify one or more sequences to include Insertions or Deletions
#'
#' @details
#' This is a lower-level function relied on by both [transmogrify()] and
#' [genomogrify()].
#'
#' Takes an [Biostrings::XString] or [Biostrings::XStringSet] object and
#' modifies the sequence to incorporate InDels.
#' The expected types of data determine the behaviour, with the following
#' expectations describing how the function will incorporate data
#'
#' | Input Data Type | Exons Required | Use Case | Returned |
#' | --------------- | -------------- | -------- | -------- |
#' | XString         | Y | Modify a Reference Transcriptome | XString    |
#' | DNAStringSet    | N | Modify a Reference Genome | DNAStringSet |
#' | BSgenome        | N | Modify a Reference Genome | DNAStringSet |
#'
#' @return A DNAStringSet or XString object (See Details)
#'
#' @param x Sequence of class XString
#' @param exons GRanges object containing exon structure for `x`
#' @param indels GRanges object with InDel locations and the alternate allele
#' @param alt_col Column containing the alternate allele
#' @param mc.cores Number of cores to use when calling [parallel::mclapply]
#' internally
#' @param verbose logical(1) Print all messages
#' @param names passed to [BSgenome::getSeq] when x is a BSgenome object
#' @param ... Passed to [parallel::mclapply]
#'
#' @seealso [transmogrify()] [genomogrify()]
#'
#' @examples
#' ## Start with a DNAStringSet
#' library(GenomicRanges)
#' seq <- DNAStringSet(c(seq1 = "AATCTGCGC"))
#' ## Define an Insertion
#' var <- GRanges("seq1:1")
#' var$ALT <- "AAA"
#' seq
#' indelcator(seq, var)
#'
#' ## To modify a single transcript
#' library(GenomicFeatures)
#' ex <- GRanges(c("seq1:1-3:+", "seq1:7-9:+"))
#' orig <- extractTranscriptSeqs(seq, GRangesList(tx1 = ex))[["tx1"]]
#' orig
#' indelcator(orig, var, exons = ex)
#'
#' @export
#' @name indelcator
#' @rdname indelcator-methods
setGeneric("indelcator", function(x, indels, ...) standardGeneric("indelcator"))
#'
#'
#' @import Biostrings
#' @import GenomicRanges
#' @importFrom methods is as
#' @importFrom IRanges subsetByOverlaps findOverlaps
#' @importFrom S4Vectors mcols mcols<- queryHits subjectHits DataFrame
#' @importFrom stats aggregate
#' @rdname indelcator-methods
#' @aliases indelcator
#' @export
setMethod(
    "indelcator",
    signature = signature(x = "XString", indels = "GRanges"),
    function(x, indels, exons, alt_col = "ALT", ...) {

        ## Get the variants matching the exons
        indels <- subsetByOverlaps(indels, exons)
        strand(indels) <- "*"
        if (length(indels) == 0) return(x)

        ## Check the sequence matches the exons
        n <- length(x)
        stopifnot(n == sum(width(exons)))

        ## Add the ID column to the variants & check the alt column
        indels <- .checkAlts(indels, alt_col)
        indels$ID <- paste0("V", seq_along(indels))
        mcols(indels) <- mcols(indels)[c("ID", alt_col)]

        ## Now set the sequence order & strandedness
        neg_stranded <- any(strand(exons) == "-")
        i <- seq_len(n)
        if (neg_stranded) i <- rev(i)

        ## Form a GPos object with position and overlapping variants
        gpos <- GPos(sort(exons, ignore.strand = TRUE))
        mcols(gpos) <- DataFrame(i = i, ID = "")
        ol <- findOverlaps(gpos, indels)
        gpos$ID[queryHits(ol)] <- indels$ID[subjectHits(ol)]
        ## Figure out where a variant starts & ends
        gpos$group <- cumsum(c(FALSE, gpos$ID[-1] != gpos$ID[-n]))

        ## Setup as a df then place each run of 'i' as a list column
        ## Testing using Views takes about 130% of the time
        df <- as.data.frame(gpos)[c("i", "ID", "group")]
        df <- df[order(df$i),]
        # df <- chop(df, i) # Is there a base version of this?
        df$ID <- factor(df$ID, levels = unique(df$ID))
        df$group <- factor(df$group, levels = unique(df$group))
        fm <- as.formula("i ~ ID + group")
        df <- aggregate(fm, data = df, FUN = list)
        df$ID <- as.character(df$ID)
        df$group <- as.integer(as.character(df$group))
        ## Start with the ref alleles
        df$seq <- lapply(df$i, \(i) as.character(x[i]))
        ## Now define the alternates
        df$alt <- df$seq
        alt <- mcols(indels)[[alt_col]]
        names(alt) <- indels$ID
        cl <- class((x)) # Needed for coercion at the final step
        if (neg_stranded) {
            new_cl <- paste0(cl, "Set")
            alt <- as.character(reverseComplement(as(alt, new_cl)))
        }

        stopifnot(all(sort(setdiff(df$ID, "") ) == sort(names(alt))))
        rows2sub <- df$ID != ""
        df$alt[rows2sub] <- alt[df$ID[rows2sub]]
        new_seq <- as(paste(df$alt, collapse = ""), cl)
        new_seq

    }
)
#'
#' @import Biostrings
#' @import GenomicRanges
#' @importFrom GenomeInfoDb seqnames seqinfo seqlevelsInUse
#' @importFrom S4Vectors splitAsList mcols mcols<-
#' @importFrom IRanges width Views start end width<-
#' @importFrom parallel mclapply
#' @rdname indelcator-methods
#' @aliases indelcator
#' @export
setMethod(
    "indelcator",
    signature(x = "DNAStringSet", indels = "GRanges"),
    function(x, indels, alt_col = "ALT", mc.cores = 1, verbose = TRUE, ...) {

        sq <- seqinfo(x)
        seq2_mod <- unique(as.character(seqnames(indels)))
        seq2_mod <- intersect(seq2_mod, seqnames(sq))
        if (length(seq2_mod) == 0) {
            message("No InDels match the reference DNAStringSet")
            return(x)
        }
        gr <- subset(GRanges(sq), seqnames %in% seq2_mod)
        grl <- splitAsList(gr, seqlevelsInUse(gr))

        indels <- .checkAlts(indels, alt_col)
        indels$deletion <- width(indels) > nchar(mcols(indels)[[alt_col]])
        indels$insertion <- width(indels) < nchar(mcols(indels)[[alt_col]])
        stopifnot(all(width(indels)[indels$insertion] == 1))
        stopifnot(all(width(indels)[indels$deletion] > 1))
        indels <- subset(indels, deletion | insertion) # remove any SNPs
        mcols(indels)$indel <- TRUE
        mcols(indels) <- mcols(indels)[c(alt_col, "indel")]
        indels <- splitAsList(indels, f = as.character(seqnames(indels)))
        x[seq2_mod] <- mclapply(
            seq2_mod,
            function(i){
                if (verbose) message(
                    "Updating ", i, "; Original length: ", length(x[[i]]),
                    appendLF = FALSE
                )
                unch <- setdiff(grl[[i]], indels[[i]])
                unch$indel <- FALSE
                all_rng <- sort(c(indels[[i]], unch))
                seq_views <- Views(x[[i]], start(all_rng), end(all_rng))
                ## Indels
                alts <- mcols(subset(all_rng, indel))[[alt_col]]
                width(seq_views)[all_rng$indel] <- nchar(alts)
                new_seq <- vector("list", length(all_rng))
                not_indel <- !all_rng$indel
                new_seq[not_indel] <- as.list(DNAStringSet(seq_views)[not_indel])
                new_seq[!not_indel] <- as.list(DNAStringSet(alts))
                out <- unlist(DNAStringSet(new_seq))
                if (verbose) message("; Updated length: ", length(out))
                out
            }, mc.cores = mc.cores, ...
        )
        x
    }
)
#' @import Biostrings
#' @import GenomicRanges
#' @importClassesFrom BSgenome BSgenome
#' @rdname indelcator-methods
#' @aliases indelcator
#' @export
setMethod(
    "indelcator",
    signature(x = "BSgenome", indels = "GRanges"),
    function(x, indels, alt_col = "ALT", mc.cores = 1, names, ...) {
        seq <- as(getSeq(x, names), "DNAStringSet")
        if (!missing(names)) names(seq) <- names
        indelcator(seq, indels, alt_col, mc.cores, ...)
    }
)
