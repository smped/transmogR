#' @title Mogrify a genome using a set of variants
#' @description
#' Use a set of insertions and deletions to modify a reference genome
#'
#' @param x A DNAStringSet or BSgenome
#' @param var GRanges object containing the variants
#' @param alt_col The name of the column with `var` containing alternate bases
#' @param seq Optional subset of sequences to apply the function to
#' @param tag Optional tag to add to all sequence names
#' @param sep Separator to place between seqnames names & tag
#' @param ... Passed to armIndello
#'
#' @return XStringSet with variant modified sequences
#'
#' @examples
#' library(GenomicRanges)
#' seq <- DNAStringSet(c(chr1 = "ACGT"))
#' var <- GRanges(c("chr1:1", "chr1:3"))
#' var$ALT <- c("C", "GG")
#' seq
#' genoMog(seq, var, tag = "mod")
#'
#'
#' @export
#' @name genoMog
#' @rdname genoMog-methods
setGeneric(
    "genoMog", function(x, var, ...){standardGeneric("genoMog")}
)
#' @importFrom S4Vectors mcols
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom IRanges width
#' @export
#' @rdname genoMog-methods
#' @aliases genoMog-methods
setMethod(
    "genoMog",
    signature = signature(x = "XStringSet", var = "GRanges"),
    function(x, var, alt_col = "ALT", seq, tag = NULL, sep = "_", ...) {

        ## 1. Identify SNPs within 'var'
        ## 2. Idenfiy InDels within 'var'
        ## 3. Substitute the SNPs
        ## 4. Substitute the InDels
        ## 5. Optionally tag sequence names
        ##    - Use a common tag + specific tags for SNPs/Insertions/Deletions

        ## Checks, then subset x
        alt_col <- match.arg(alt_col, colnames(mcols(var)))
        if (missing(seq)) seq <- seqlevels(x)
        if (any(duplicated(seq))) stop("Duplicate sequence names not permitted")
        seq <- intersect(seq, seqlevels(x))
        if (is.null(seq)) return(x)
        stopifnot(any(seq %in% seqlevels(var)))
        x <- x[seq]

        ## Deal with any NA values
        is_na <- is.na(mcols(var)[[alt_col]])
        if (any(is_na)) message("NA values found in alt_col. Skipping these.")
        var <- var[!is_na]
        ## Separate into snps & indels
        is_snp <- width(var) == 1 & nchar(mcols(var)[[alt_col]]) == 1
        snps <- var[is_snp]
        indels <- var[!is_snp]

        ## Overwrite SNPs
        new_seq <- owl(x, snps, alt_col)
        ## Arm with Insertions/Deletions
        new_seq <- armIndello(new_seq, indels, ...)

        ## Sort out tags for sequence names
        new_names <- names(new_seq)
        if (!is.null(tag)) new_names <- paste(new_names, tag, sep = sep)
        names(new_seq) <- new_names
        new_seq

    }
)
#' @importFrom Biostrings getSeq
#' @export
#' @rdname genoMog-methods
#' @aliases genoMog-methods
setMethod(
    "genoMog",
    signature = signature(x = "BSgenome", var = "GRanges"),
    function(x, var, alt_col = "ALT", seq, tag = NULL, sep = "_", ...) {
        ## Setup the sequence info
        message("Extracting sequences as a DNAStringSet...", appendLF = FALSE)
        x <- getSeq(x)
        message("done")
        genoMog(x, var, alt_col, seq, tag, sep, ...)
    }
)
