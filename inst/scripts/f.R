f <- function(x, indels, alt_col = "ALT", mc.cores = 1, verbose = TRUE, ...) {

    sq <- seqinfo(x)
    seq2_mod <- unique(as.character(seqnames(indels)))
    seq2_mod <- intersect(seq2_mod, seqnames(sq))
    if (length(seq2_mod) == 0) {
        message("No InDels match the reference DNAStringSet")
        return(x)
    }
    gr <- subset(GRanges(sq), seqnames %in% seq2_mod)
    grl <- splitAsList(gr, seqlevelsInUse(gr))
    cl <- as.character(class(x))

    browser()

    indels <- splitAsList(indels, f = as.character(seqnames(indels)))

    # x[seq2_mod] <-
    mclapply(
        seq2_mod,
        \(i) {
            if (verbose) message(
                "Updating ", i, "; Original length: ", length(x[[i]]),
                appendLF = FALSE
            )
            unch <- setdiff(grl[[i]], indels[[i]])
            all_rng <- sort(c(indels[[i]], unch))
            ol <- findOverlaps(all_rng, indels[[i]])
            seq_views <- Views(x[[i]], start(all_rng), end(all_rng))
            new_seq <- as(seq_views, cl)

        }
    )



    indels <- .checkAlts(indels, alt_col) # What's this doing?
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

            new_seq <- DNAStringSet(seq_views)
            not_indel <- !all_rng$indel
            new_seq[!not_indel] <- DNAStringSet(alts)
            out <- unlist(DNAStringSet(new_seq))
            if (verbose) message("; Updated length: ", length(out))
            out

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
