#' @title Mogrify a genome using a set of variants
#' @description
#' Use a set of SNPS, insertions and deletions to modify a reference genome
#'
#' @param x A DNAStringSet or BSgenome
#' @param var GRanges object containing the variants
#' @param alt_col The name of the column with `var` containing alternate bases
#' @param names Sequence names to be mogrified
#' @param tag Optional tag to add to all sequence names which were modified
#' @param sep Separator to place between seqnames names & tag
#' @param var_tags logical(1) Add tags indicating which type of variant were
#' incorporated, with 's', 'i' and 'd' represeenting SNPs, Insertions and
#' Deletions respectively
#' @param var_sep Separator between any previous tags and variant tags
#' @param ... Passed to armIndello
#'
#' @return XStringSet with variant modified sequences
#'
#' @examples
#' library(GenomicRanges)
#' dna <- DNAStringSet(c(chr1 = "ACGT", chr2 = "AATTT"))
#' var <- GRanges(c("chr1:1", "chr1:3", "chr2:1-3"))
#' var$ALT <- c("C", "GG", "A")
#' dna
#' mogrifyGenome(dna, var)
#' mogrifyGenome(dna, var, tag = "mod")
#' mogrifyGenome(dna, var, var_tags = TRUE)
#'
#'
#' @export
#' @name mogrifyGenome
#' @rdname mogrifyGenome-methods
setGeneric(
    "mogrifyGenome", function(x, var, ...){standardGeneric("mogrifyGenome")}
)
#' @importFrom S4Vectors mcols
#' @importFrom GenomeInfoDb seqlevels seqnames
#' @importFrom IRanges width
#' @export
#' @rdname mogrifyGenome-methods
#' @aliases mogrifyGenome-methods
setMethod(
    "mogrifyGenome",
    signature = signature(x = "XStringSet", var = "GRanges"),
    function(
        x, var, alt_col = "ALT", tag = NULL, sep = "_",
        var_tags = FALSE, var_sep = "_", ...
    ) {

        ## 1. Identify SNPs within 'var'
        ## 2. Idenfiy InDels within 'var'
        ## 3. Substitute the SNPs
        ## 4. Substitute the InDels
        ## 5. Optionally tag sequence names
        ##    - Use a common tag + specific tags for SNPs/Insertions/Deletions

        ## Check the variants are valid
        var <- .checkAlts(var, alt_col)
        ## Separate into snps & indels
        var <- subset(var, seqnames %in% seqlevels(x))
        if (length(var) == 0) return(x)
        is_snp <- width(var) == 1 & nchar(mcols(var)[[alt_col]]) == 1
        snps <- var[is_snp]
        indels <- var[!is_snp]

        ## Overwrite SNPs
        new_seq <- owl(x, snps, alt_col)
        ## Arm with Insertions/Deletions
        new_seq <- armIndello(new_seq, indels, ...)

        ## Sort out tags for sequence names which were modified
        if (!is.null(tag)) {
            ol <- names(new_seq) %in% seqlevels(var)
            names(new_seq)[ol] <- paste(names(new_seq)[ol], tag, sep = sep)
        }
        if (var_tags) {
            ol_snp <- names(x) %in% seqnames(snps)
            ol_ins <- names(x) %in% seqnames(subset(indels, width == 1))
            ol_del <- names(x) %in% seqnames(subset(indels, width > 1))
            any_var <- ol_snp | ol_ins | ol_del
            suff <- rep_len("", length(new_seq))
            suff[any_var] <- var_sep
            suff[ol_snp] <- paste0(suff[ol_snp], "s")
            suff[ol_ins] <- paste0(suff[ol_ins], "i")
            suff[ol_del] <- paste0(suff[ol_del], "d")
            names(new_seq) <- paste0(names(new_seq), suff)
        }
        new_seq

    }
)
#' @importFrom Biostrings getSeq
#' @export
#' @rdname mogrifyGenome-methods
#' @aliases mogrifyGenome-methods
setMethod(
    "mogrifyGenome",
    signature = signature(x = "BSgenome", var = "GRanges"),
    function(
        x, var, alt_col = "ALT", names, tag = NULL, sep = "_",
        var_tags = FALSE, var_sep = "_", ...
    ) {
        ## Setup the sequence info
        message("Extracting sequences as a DNAStringSet...", appendLF = FALSE)
        x <- getSeq(x, names)
        message("done")
        mogrifyGenome(x, var, alt_col, tag, sep, var_tags = FALSE, var_sep = "_", ...)
    }
)