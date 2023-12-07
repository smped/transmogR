#' @title Mogrify a transcriptome using a set of variants
#' @description
#' Use a set of SNPs, insertions and deletions to modify a reference
#' transcriptome
#'
#' @details
#' Produce a set of variant modified transcript sequences from a standard
#' reference genome.
#' Supported variants are SNPs, Insertions and Deletions
#'
#' Ranges needing to be masked, such as the Y-chromosome, or Y-PAR can be
#' provided.
#'
#' **It should be noted that this is a time consuming process**
#' Inclusion of a large set of insertions and deletions across an entire
#' transcriptome can involve individually modifying many thousands of
#' transcripts, which can be a computationally demanding task.
#' Whilst this can be parallelised using an appropriate number of cores, this
#' may also prove taxing for lower power laptops, and pre-emptively closing
#' memory hungry programs such as Slack, or internet browers may be prudent.
#'
#' @param x Reference genome as either a DNAStringSet or BSgenome
#' @param var GRanges object containing the variants
#' @param exons GRanges object with ranges representing exons
#' @param alt_col Column from `var` containing alternate bases
#' @param trans_col Column from 'exons' containing the transcript_id
#' @param omit_ranges GRanges object containing ranges to omit, such as PAR-Y
#' regions, for example
#' @param tag Optional tag to add to all sequence names which were modified
#' @param sep Separator to place between seqnames names & tag
#' @param var_tags logical(1) Add tags indicating which type of variant were
#' incorporated, with 's', 'i' and 'd' representing SNPs, Insertions and
#' Deletions respectively
#' @param var_sep Separator between any previous tags and variant tags
#' @param verbose logical(1) Include informative messages, or operate silently
#' @param mc_cores Number of cores to be used when multi-threading via
#' [parallel::mclapply]
#' @param ... Not used
#'
#' @return An XStringSet
#'
#' @examples
#' library(GenomicRanges)
#' library(GenomicFeatures)
#' seq <- DNAStringSet(c(chr1 = "ACGTAAATGG"))
#' exons <- GRanges(c("chr1:1-3:-", "chr1:7-9:-"))
#' exons$transcript_id <- c("trans1")
#'
#' # When using extractTranscriptSeqs -stranded exons need to be sorted by end
#' exons <- sort(exons, decreasing = TRUE, by = ~end)
#' exons
#' trByExon <- splitAsList(exons, exons$transcript_id)
#'
#' # Check the sequences
#' seq
#' extractTranscriptSeqs(seq, trByExon)
#'
#' # Define some variants
#' var <- GRanges(c("chr1:2", "chr1:8"))
#' var$ALT <- c("A", "GGG")
#'
#' # Include the variants adding tags to indicate a SNP and indel
#' # The exons GRanges object will be split by transcript internally
#' mogrifyTranscriptome(seq, var, exons, var_tags = TRUE)
#'
#'
#' @export
#' @name mogrifyTranscriptome
#' @rdname mogrifyTranscriptome-methods
setGeneric(
    "mogrifyTranscriptome",
    function(x, var, exons, ...){standardGeneric("mogrifyTranscriptome")}
)
#' @importFrom S4Vectors mcols splitAsList
#' @importFrom GenomeInfoDb seqlevels seqnames
#' @importFrom IRanges width subsetByOverlaps
#' @importFrom GenomicFeatures extractTranscriptSeqs
#' @importFrom parallel mclapply
#' @importFrom stats as.formula
#' @export
#' @rdname mogrifyTranscriptome-methods
#' @aliases mogrifyTranscriptome-methods
setMethod(
    "mogrifyTranscriptome",
    signature = signature(x = "XStringSet", var = "GRanges", exons = "GRanges"),
    function(
        x, var, exons, alt_col = "ALT", trans_col = "transcript_id",
        omit_ranges = NULL, tag = NULL, sep = "_",
        var_tags = FALSE, var_sep = "_", verbose = TRUE, mc_cores = 1, ...
    ) {

        ## 1. Identify SNPs within 'var'
        ## 2. Idenfiy InDels within 'var'
        ## 3. Substitute the SNPs into the reference
        ## 4. Extract Transcripts
        ## 5. Substitute the InDels
        ## 6. Optionally tag sequence names
        ##    - Use a common tag + specific tags for SNPs/Insertions/Deletions

        ## Checks
        trans_col <- match.arg(trans_col, colnames(mcols(exons)))
        var <- .checkAlts(var, alt_col)

        ## Separate into snps & indels
        var <- subset(var, seqnames %in% seqlevels(x))
        if (length(var) == 0) return(x)
        is_snp <- width(var) == 1 & nchar(mcols(var)[[alt_col]]) == 1
        snps <- var[is_snp]
        indels <- var[!is_snp]

        ## Remove any unwanted transcripts
        if (!is.null(omit_ranges)) {
            stopifnot(is(omit_ranges, "GRanges"))
            omit <- mcols(subsetByOverlaps(exons, omit_ranges))[[trans_col]]
            exons <- exons[!mcols(exons)[[trans_col]] %in% omit]
        }
        ## Sort by strand. Exons cannot be unstranded
        strand <- strand(exons)
        if (any(strand == "*"))
            stop("Unstranded exons found. Exons must be stranded")
        exList <- splitAsList(exons, strand)
        exList[["+"]] <- sort(exList[["+"]])
        fm <- as.formula("~seqnames + end")
        exList[["-"]] <- sort(exList[["-"]], decreasing = TRUE, by = fm)
        exons <- unlist(exList)

        ## Find those with InDels as we only need to really work on these
        trans_with_any <- mcols(subsetByOverlaps(exons, var))[[trans_col]]
        trans_with_indel <- mcols(subsetByOverlaps(exons, indels))[[trans_col]]
        trans_with_indel <- unique(trans_with_indel)
        if (verbose) message(
            length(trans_with_indel), " transcripts found with indels"
        )

        ## Modify the genome & extract sequences including SNPs
        new_ref <- owl(x, snps, alt_col = alt_col)
        ex_by_trans <- splitAsList(exons, mcols(exons)[[trans_col]])
        all_seq <- extractTranscriptSeqs(new_ref, ex_by_trans)

        ## Now modify the transcripts with InDels
        cl <- class(all_seq)
        new_trans_seq <- mclapply(
            trans_with_indel,
            function(id) {
                armIndello(
                    all_seq[[id]], indels, ex_by_trans[[id]], alt_col
                )
            }, mc.cores = mc_cores
        )
        names(new_trans_seq) <- trans_with_indel
        new_trans_seq <- as(new_trans_seq, cl)
        all_seq[trans_with_indel] <- new_trans_seq[trans_with_indel]

        ## Add tags where needed
        if (!is.null(tag)) {
            ol <- names(all_seq) %in% trans_with_any
            names(all_seq)[ol] <- paste(names(all_seq)[ol], tag, sep = sep)
        }
        if (var_tags) {
            suff <- rep_len("", length(ex_by_trans))
            trans_with_snp <- names(subsetByOverlaps(ex_by_trans, snps))
            trans_with_ins <- names(
                subsetByOverlaps(ex_by_trans, subset(indels, width == 1))
            )
            trans_with_del <- names(
                subsetByOverlaps(ex_by_trans, subset(indels, width > 1))
            )
            suff[names(ex_by_trans) %in% trans_with_any] <- var_sep
            suff[names(ex_by_trans) %in% trans_with_snp] <- paste0(suff, "s")
            suff[names(ex_by_trans) %in% trans_with_ins] <- paste0(suff, "i")
            suff[names(ex_by_trans) %in% trans_with_del] <- paste0(suff, "d")
            names(all_seq) <- paste0(names(all_seq), suff)
        }
        all_seq

    }
)
#' @importFrom Biostrings getSeq
#' @importFrom GenomeInfoDb seqnames
#' @importFrom GenomicRanges granges
#' @export
#' @rdname mogrifyTranscriptome-methods
#' @aliases mogrifyTranscriptome-methods
setMethod(
    "mogrifyTranscriptome",
    signature = signature(x = "BSgenome", var = "GRanges", exons = "GRanges"),
    function(
        x, var, exons, alt_col = "ALT", trans_col = "transcript_id",
        omit_ranges = NULL, tag = NULL, sep = "_",
        var_tags = FALSE, var_sep = "_", verbose = TRUE, mc_cores = 1, ...
    ) {
        ## Setup the sequence info, only extracting those with a transcript
        all_gr <- sort(c(granges(var), granges(exons)))
        seq_to_get <- unique(seqnames(all_gr))
        if (verbose) message(
            "Extracting ", length(seq_to_get),
            " sequences as a DNAStringSet...", appendLF = FALSE
        )
        x <- as(getSeq(x, seq_to_get), "DNAStringSet")
        if (verbose) message("done")
        mogrifyTranscriptome(
            x, var, exons, alt_col, trans_col, omit_ranges, tag, sep, var_tags,
            var_sep, verbose, mc_cores, ...
        )
    }
)
