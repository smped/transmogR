#' @title Substitute InDels into one or more sequences
#'
#' @description
#' Modify one or more sequences to include Insertions or Deletions
#'
#' @details
#' Takes an [Biostrings::XString] or [Biostrings::XStringSet] object and moddfies
#' the sequence to incorporate InDels.
#' The expected types of data determine the behaviour, with the following
#' expectations describing how the function will incorporate data
#'
#' | Input Data Type | Exons Required | Use Case | Returned |
#' | --------------- | -------------- | -------- | -------- |
#' | XString         | Y | Modify a Reference Transcriptome | XString    |
#' | XStringSet      | Y | Modify a Reference Transcriptome | XStringSet |
#' | DNAStringSet    | N | Modify a Reference Genome | DNAStringSet |
#' | BSgenome        | N | Modify a Reference Genome | DNAStringSet |
#'
#'
#' @param x Sequence of class XString
#' @param exons GRanges object containing exon structure for `x`
#' @param indels GRanges object with InDel locations and the alternate allele
#' @param alt_col Column containing the alternate allele
#' @param mc_cores Number of cores to use when calling [parallel::mclapply]
#' internally
#' @param names passed to [BSgenome::getSeq] when x is a BSgenome object
#' @param ... Not used
#'
#'
#' @export
#' @name armIndello
#' @rdname armIndello-methods
setGeneric(
  "armIndello", function(x, indels, exons, ...){standardGeneric("armIndello")}
  # Originally subInDel
)
#'
#'
#' @import Biostrings
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom methods is as
#' @importFrom IRanges subsetByOverlaps findOverlaps
#' @importFrom S4Vectors mcols 'mcols<-' queryHits subjectHits DataFrame
#' @importFrom GenomicRanges strand 'strand<-' GPos
#' @importFrom tidyr chop
#' @rdname armIndello-methods
#' @aliases armIndello
#' @export
setMethod(
  "armIndello",
  signature = signature(x = "XString", indels = "GRanges", exons = "GRanges"),
  function(x, indels, exons, alt_col = "ALT", ...) {

    ## Check classes (Not needed for S4)
    cl <- class((x))

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

    ## Form a GPos object with positional information and overlapping variants
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
    df <- chop(df, i)
    ## Start with the ref alleles
    df$seq <- lapply(df$i, \(i) as.character(x[i]))
    ## Now define the alternates
    df$alt <- df$seq
    alt <- mcols(indels)[[alt_col]]
    names(alt) <- indels$ID
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
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom methods as
#' @rdname armIndello-methods
#' @aliases armIndello
#' @export
setMethod(
  "armIndello",
  signature = signature(x = "XStringSet", indels = "GRanges", exons = "GRanges"),
  function(x, indels, exons, alt_col = "ALT", ...) {
    cmn_seq <- intersect(seqnames(x), seqnames(indels))
    cmn_seq <- intersect(cmn_seq, seqnames(exons))
    if (length(cmn_seq) == 0) stop("No shared sequences found")
    cl <- class(x)
    x <- unlist(x) ## Coerce to an XString object
    out <- armIndello(x, indels, exons, alt_col, ...)
    as(out, cl)
  }
)
#'
#' @import Biostrings
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb seqnames seqinfo seqlevelsInUse
#' @importFrom S4Vectors splitAsList mcols 'mcols<-'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges width Views start end 'width<-'
#' @importFrom parallel mclapply
#' @rdname armIndello-methods
#' @aliases armIndello
#' @export
setMethod(
  "armIndello",
  signature(x = "DNAStringSet", indels = "GRanges", exons = "missing"),
  function(x, indels, exons, alt_col = "ALT", mc_cores = 1, ...) {
    ## Still can't run an entire genome on the laptop in parallel.
    ## Takes 3 mins using SerialParam() which is pretty good
    # armIndello(hg38_mod[1:2], indel = gr_indel, BPPARAM = MulticoreParam(2))
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
    mcols(indels) <- mcols(indels)[c(alt_col, "deletion", "insertion")]
    indels <- splitAsList(indels, f = as.character(seqnames(indels)))
    x[seq2_mod] <- mclapply(
      seq2_mod,
      function(i){
        message(
          "Updating ", i, "; Original length: ", length(x[[i]]),
          appendLF = FALSE
        )
        unch <- setdiff(grl[[i]], indels[[i]])
        unch$deletion <- FALSE
        unch$insertion <- FALSE
        all_rng <- sort(c(indels[[i]], unch))
        seq_views <- Views(x[[i]], start(all_rng), end(all_rng))
        ## Deletions
        width(seq_views)[all_rng$deletion] <- 1
        ## Insertions
        alts <- mcols(subset(all_rng, insertion))[[alt_col]]
        width(seq_views)[all_rng$insertion] <- nchar(alts)
        new_seq <- vector("list", length(all_rng))
        not_ins <- !all_rng$insertion
        new_seq[not_ins] <- as.list(DNAStringSet(seq_views)[not_ins])
        new_seq[!not_ins] <- as.list(DNAStringSet(alts))
        out <- unlist(DNAStringSet(new_seq))
        message("; Updated length: ", length(out))
        out
      }, mc.cores = mc_cores
    )
    x
  }
)
#' @import Biostrings
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom BSgenome BSgenome
#' @rdname armIndello-methods
#' @aliases armIndello
#' @export
setMethod(
  "armIndello",
  signature(x = "BSgenome", indels = "GRanges", exons = "missing"),
  function(x, indels, exons, alt_col = "ALT", mc_cores = 1, names, ...) {
    seq <- getSeq(x, names)
    armIndello(seq, indels, exons, alt_col, mc_cores, ...)
  }
)
