#' @title Show Variants by Impacted Columns
#'
#' @description
#' Produce an UpSet plot showing unique values from a given column
#'
#' @details
#' Take a set of variants, classify them as SNV, Insertion and Deletion, then
#' using a GRanges object, produce an UpSet plot showing impacted values from
#' a given column
#'
#' @return An UpSet plot
#'
#' @seealso [ComplexUpset::upset]
#'
#' @param gr GRanges object with ranges representing a key feature such as exons
#' @param var GRanges object with variants in a given column
#' @param alt_col Column within `var` containing the alternate allele
#' @param mcol The column within `gr` to summarise results by
#' @param ... Passed to [ComplexUpset::upset]
#' @param intersection_args See [ComplexUpset::intersection_size] for possible
#' values
#' @param intersection_lab Y-axis label for the intersection panel
#' @param set_geom Passed to [ComplexUpset::upset_set_size]
#' @param set_expand Expand the set-size axis by this amount
#' @param set_counts logical(1) Show counts on set sizes
#' @param hjust_counts Horizontal adjustment of counts, if being shown
#' @param set_lab X-axis label for the set-sizes panel
#' @param title Summary title to show above the intersection panel. Can be
#' hidden by setting to NULL
#'
#' @examples
#' library(rtracklayer)
#' library(VariantAnnotation)
#' gtf <- import.gff(
#'     system.file("extdata/gencode.v44.subset.gtf.gz", package = "transmogR"),
#'     feature.type = "exon"
#' )
#' vcf <- system.file("extdata/1000GP_subset.vcf.gz", package = "transmogR")
#' var <- rowRanges(readVcf(vcf, param = ScanVcfParam(fixed = "ALT")))
#' upsetVarByCol(gtf, var)
#'
#'
#'
#' @importFrom S4Vectors splitAsList mcols
#' @importFrom IRanges subsetByOverlaps
#' @importFrom ComplexUpset intersection_size upset_set_size upset
#' @importFrom ggplot2 aes geom_bar geom_text after_stat ggtitle position_stack
#' @importFrom ggplot2 scale_y_reverse scale_y_continuous expansion ylab
#' @importFrom rlang list2 ':='
#' @export
upsetVarByCol <- function(
        gr, var, alt_col = "ALT", mcol = "transcript_id", ...,
        intersection_args = list(), intersection_lab = "Intersection Size",
        set_geom = geom_bar(width = 0.6), set_expand = 0.2, set_counts = TRUE,
        hjust_counts = 1.1, set_lab = "Set Size", title
){

    stopifnot(is(gr, "GRanges"))
    mcol <- match.arg(mcol, colnames(mcols(gr)))

    stopifnot(is(var, "GRanges"))
    alt_col <- match.arg(alt_col, colnames(mcols(var)))
    var$type <- calvInDel(var, alt_col)
    var_list <- splitAsList(var, var$type)

    ids <- unique(mcols(gr)[[mcol]])
    ol <- lapply(var_list, \(x) mcols(subsetByOverlaps(gr, x))[[mcol]])
    ol <- lapply(ol, unique)
    ol_ids <- unique(unlist(ol))
    set_max <- max(vapply(ol, length, integer(1)))

    if (missing(title)) title <- paste0(
        scales::comma(length(ol_ids), 1), " Overlap; ",
        scales::comma(length(setdiff(ids, ol_ids)), 1),
        " No Overlap (", mcol, ")"
    )

    df_list <- list(ids = ol_ids, lapply(ol, \(x) ol_ids %in% x))
    df <- as.data.frame(df_list)
    base_args <- .makeIntersectionArgs(intersection_args)
    base_ann <- do.call("intersection_size", base_args) +
        scale_y_continuous(expand = expansion(c(0, 0.05)))
    if (is.character(title)) base_ann <- base_ann + ggtitle(title)

    sets <- upset_set_size(geom = set_geom) + ylab(set_lab) +
        scale_y_reverse(expand = expansion(c(set_expand, 0)))
    if (set_counts) {
        count <- c()
        sets <- sets + geom_text(
            aes(label = after_stat(count)), hjust = hjust_counts, stat = 'count'
        )
    }

    upset(
        df, names(var_list),
        base_annotations = list2("{intersection_lab}" := base_ann),
        set_sizes = sets, ...
    )

}


