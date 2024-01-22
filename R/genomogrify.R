#' @title Mogrify a genome using a set of variants
#'
#' @description
#' Use a set of SNPS, insertions and deletions to modify a reference genome
#'
#' @details
#' This function is designed to create a variant-modified reference genome,
#' intended to be included as a set of decoys when using salmon in selective
#' alignment mode.
#' Sequence lengths will change if InDels are included and any coordinate-based
#' information will be lost on the output of this function.
#'
#' Tags are able to be added to any modified sequence to assist identifying any
#' changes that have been made to a sequence.
#'
#'
#' @param x A DNAStringSet or BSgenome
#' @param var GRanges object containing the variants, or a
#' [VariantAnnotation::VcfFile]
#' @param alt_col The name of the column with `var` containing alternate bases
#' @param mask Optional GRanges object defining regions to be masked with an 'N'
#' @param names Sequence names to be mogrified
#' @param tag Optional tag to add to all sequence names which were modified
#' @param sep Separator to place between seqnames names & tag
#' @param var_tags logical(1) Add tags indicating which type of variant were
#' incorporated, with 's', 'i' and 'd' representing SNPs, Insertions and
#' Deletions respectively
#' @param var_sep Separator between any previous tags and variant tags
#' @param which GRanges object passed to [VariantAnnotation::ScanVcfParam] if
#' using a VCF directly
#' @param verbose logical(1) Print progress messages while running
#' @param ... Passed to [parallel::mclapply]
#'
#' @return XStringSet with variant modified sequences
#'
#' @examples
#' library(GenomicRanges)
#' dna <- DNAStringSet(c(chr1 = "ACGT", chr2 = "AATTT"))
#' var <- GRanges(c("chr1:1", "chr1:3", "chr2:1-3"))
#' var$ALT <- c("C", "GG", "A")
#' dna
#' genomogrify(dna, var)
#' genomogrify(dna, var, tag = "mod")
#' genomogrify(dna, var, var_tags = TRUE)
#' genomogrify(dna, var, mask = GRanges("chr2:1-5"), var_tags = TRUE)
#'
#'
#' @export
#' @name genomogrify
#' @rdname genomogrify-methods
setGeneric(
    "genomogrify", function(x, var, ...){standardGeneric("genomogrify")}
)
#' @importFrom S4Vectors mcols
#' @importFrom GenomeInfoDb seqlevels seqnames
#' @importFrom IRanges width
#' @export
#' @rdname genomogrify-methods
#' @aliases genomogrify-methods
setMethod(
    "genomogrify",
    signature = signature(x = "XStringSet", var = "GRanges"),
    function(
        x, var, alt_col = "ALT", mask, tag = NULL, sep = "_",
        var_tags = FALSE, var_sep = "_", verbose = TRUE, ...
    ) {

        ## 1. Identify SNPs within 'var'
        ## 2. Idenfiy InDels within 'var'
        ## 3. Substitute the SNPs
        ## 4. Substitute the InDels
        ## 5. Optionally tag sequence names
        ##    - Use a common tag + specific tags for SNPs/Insertions/Deletions

        if (missing(mask)) mask <- GRanges()
        stopifnot(is(mask, "GRanges"))

        ## Check the variants are valid
        var <- .checkAlts(var, alt_col)
        var <- var[!overlapsAny(var, mask)]
        ## Separate into snps & indels
        var <- subset(var, seqnames %in% seqlevels(x))
        if (length(var) == 0) return(x)
        type <- varTypes(var, alt_col)
        snps <- var[type == "SNV"]
        indels <- var[type != "SNV"]

        ## Overwrite SNPs
        new_seq <- owl(x, snps, alt_col)
        ## Ensure masked regions are all 'N'
        mask <- subset(mask, seqnames %in% names(new_seq))
        if (length(mask) > 0) {
            if (verbose) message("Applying mask")
            mask_grl <- splitAsList(mask, as.character(seqnames(mask)))
            for (m in names(mask_grl)) {
                i <- start(GPos(mask_grl[[m]]))
                new_seq[[m]][i] <- "N"
            }
        }
        ## Add Insertions/Deletions
        new_seq <- indelcator(new_seq, indels, verbose = verbose, ...)

        ## Sort out tags for sequence names which were modified
        if (!is.null(tag)) {
            ol <- names(new_seq) %in% as.character(seqnames(var))
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
#' @rdname genomogrify-methods
#' @aliases genomogrify-methods
setMethod(
    "genomogrify",
    signature = signature(x = "BSgenome", var = "GRanges"),
    function(
        x, var, alt_col = "ALT", mask, names, tag = NULL, sep = "_",
        var_tags = FALSE, var_sep = "_", verbose = TRUE, ...
    ) {
        ## Setup the sequence info
        if (verbose) message(
            "Extracting sequences as a DNAStringSet...", appendLF = FALSE
        )
        seq <- as(getSeq(x, names), "DNAStringSet")
        if (!missing(names)) names(seq) <- names
        if (verbose) message("done")
        genomogrify(
            seq, var, alt_col, mask, tag, sep,
            var_tags = var_tags, var_sep = "_", verbose = verbose, ...
        )
    }
)
#' @importClassesFrom VariantAnnotation VcfFile
#' @export
#' @rdname genomogrify-methods
#' @aliases genomogrify-methods
setMethod(
    "genomogrify",
    signature = signature(x = "BSgenome", var = "VcfFile"),
    function(
        x, var, alt_col = "ALT", mask, names, tag = NULL, sep = "_",
        var_tags = FALSE, var_sep = "_", which, verbose = TRUE, ...
    ) {
        var <- .parseVariants(var, alt_col, which)
        genomogrify(
            x, var, alt_col, mask, names, tag, sep, var_tags = var_tags,
            var_sep = "_", verbose = verbose,  ...
        )
    }
)
#' @importClassesFrom VariantAnnotation VcfFile
#' @export
#' @rdname genomogrify-methods
#' @aliases genomogrify-methods
setMethod(
    "genomogrify",
    signature = signature(x = "XStringSet", var = "VcfFile"),
    function(
        x, var, alt_col = "ALT", mask, tag = NULL, sep = "_",
        var_tags = FALSE, var_sep = "_", which, verbose = TRUE, ...
    ) {
        var <- .parseVariants(var, alt_col, which)
        genomogrify(
            x, var, alt_col, mask, tag, sep,
            var_tags = var_tags, var_sep = "_", verbose = verbose, ...
        )
    }
)
