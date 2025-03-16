#' @title Create a set of tags indicating overlap status with variants
#'
#' @description
#' Create a set of tags indicating overlap status with variants
#'
#' @details
#' Take a GRanges or GRangesList and compare against a set of variants.
#' Variants will be classified into SNV, Insertions and Deletions using
#' [varTypes()] and tags defined.
#' An overall set of tags defining any overlap can be created by themselves.
#' An additional set of tags containing 's', 'i' or 'd' to indicate overlap
#' with an SNV, Insertion or Deletion can also be created, with the
#' concatentation of both tags being returned.
#'
#' @return Character vector of the same length as `x`
#'
#' @param x GRanges or GRangesList
#' @param var Set of variants for `x` to be compared to
#' @param tag Tag to be added for all overlapping positions
#' @param var_tags logical(1) Include 's', 'i' and 'd' tags. See details
#' @param sep Separator added between tag and var_tags
#' @param pre Separator to add at the start of returned tags
#' @param ... Passed to [cleanVariants()]
#'
#' @examples
#' # Load the included subset of 1000 Genomes Variants
#' library(VariantAnnotation)
#' vcf <- system.file("extdata/1000GP_subset.vcf.gz", package = "transmogR")
#' vcf <- VcfFile(vcf)
#' var <- cleanVariants(vcf)
#' # Now load some exons, then split by transcript, subsetting to the first 40
#' library(rtracklayer)
#' f <- system.file("extdata/gencode.v44.subset.gtf.gz", package = "transmogR")
#' gtf <- import.gff(f, feature.type = "exon")
#' exon_by_trans <- splitAsList(gtf, gtf$transcript_id)[1:40]
#' # And produce tags based on the overlapping variants within the exons
#' # Overlapping SNVs will return an 's' whilst insertions include an 'i'
#' varTags(exon_by_trans, var, tag = "1000GP")
#'
#'
#' @importFrom S4Vectors splitAsList
#' @export
varTags <- function(
        x, var, tag = NULL, var_tags = TRUE,  sep = "_", pre = sep, ...
) {

    stopifnot(is(x, "GenomicRanges_OR_GenomicRangesList"))
    ## overlapsAny behaves conveniently for both classes!
    var <- cleanVariants(var, ...)
    n <- length(x)
    # The overall tag is the prefix & the variant type tags are the suffix
    tag_prefix <- tag_suffix <- rep_len("", n)

    if (var_tags) {

        tags <- c("s", "i", "d")
        split_var <- splitAsList(var, varTypes(var))
        names(split_var) <- tolower(substr(names(split_var), 1, 1))
        tags <- intersect(tags, names(split_var)) # Avoid any missing var type
        split_var <- split_var[tags]

        ## Return a list with tags
        ol_list <- lapply(
            tags,
            \(i) {
                ret_val <- rep_len("", n)
                ret_val[overlapsAny(x, split_var[[i]])] <- i
                ret_val
            }
        )
        tag_mat <- do.call("cbind", ol_list)
        tag_suffix <- apply(tag_mat, 1, paste, collapse = "")

    }
    if (!is.null(tag)) {
        stopifnot(is.character(tag) & length(tag) == 1)
        ol <- overlapsAny(x, var)
        if (!var_tags) sep <- ""
        tag_prefix[ol] <- paste0(tag, sep)
    }
    out <- paste0(tag_prefix, tag_suffix)
    out[nchar(out) > 0] <- paste0(pre, out[nchar(out) > 0])
    out

}

