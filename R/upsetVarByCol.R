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
#' @seealso [SimpleUpset::simpleUpSet]
#'
#' @param gr GRanges object with ranges representing a key feature such as exons
#' @param var GRanges object with variants in a given column
#' @param alt_col Column within `var` containing the alternate allele
#' @param mcol The column within `gr` to summarise results by
#' @param ... Passed to [SimpleUpset::simpleUpSet]
#' @param fill Optional column in gr used to fill intersections and sets
#' @param fill_scale Discrete ggplot2 scale for filling bars. Ignored if
#' `fill = NULL`
#' @param expand_sets Expand the set-size x-axis by this amount
#' @param expand_intersect Expand the intersection y-axis by this amount
#' @param hj_sets Horizontal adjustment of set size labels
#' @param vj_intersect Vertical adjustment of intersection size labels
#' @param label_size Control the size of both intersection and sit size labels
#' @param title Summary title to show above the intersection panel. Can be
#' hidden by setting to NULL
#'
#' @examples
#' library(rtracklayer)
#' library(VariantAnnotation)
#' library(ggplot2)
#' gtf <- import.gff(
#'     system.file("extdata/gencode.v44.subset.gtf.gz", package = "transmogR"),
#'     feature.type = "exon"
#' )
#' vcf <- system.file("extdata/1000GP_subset.vcf.gz", package = "transmogR")
#' var <- rowRanges(readVcf(vcf, param = ScanVcfParam(fixed = "ALT")))
#' upsetVarByCol(gtf, var)
#' upsetVarByCol(
#'   gtf, var, fill = "transcript_type",
#'   fill_scale = scale_fill_brewer(palette = "Set1")
#' )
#'
#'
#' @importFrom S4Vectors splitAsList mcols
#' @importFrom IRanges subsetByOverlaps
#' @importFrom patchwork plot_annotation
#' @importFrom ggplot2 scale_fill_discrete guides guide_none
#' @export
upsetVarByCol <- function(
        gr, var, alt_col = "ALT", mcol = "transcript_id", ...,
        fill = NULL, fill_scale = scale_fill_discrete(),
        expand_sets = 0.2, hj_sets = 1.1,
        expand_intersect = 0.1, vj_intersect = -0.5,
        label_size = 3.5, title
){

    if (!requireNamespace('SimpleUpset', quietly = TRUE))
        stop("Please install 'SimpleUpset' to use this function.")

    stopifnot(is(gr, "GRanges"))
    mcol <- match.arg(mcol, colnames(mcols(gr)))

    stopifnot(is(var, "GRanges"))
    alt_col <- match.arg(alt_col, colnames(mcols(var)))
    var$type <- varTypes(var, alt_col)
    var_list <- splitAsList(var, var$type)

    ids <- unique(mcols(gr)[[mcol]])
    ol <- lapply(var_list, \(x) mcols(subsetByOverlaps(gr, x))[[mcol]])
    ol <- lapply(ol, unique)
    ol_ids <- unique(unlist(ol))
    if (length(ol_ids) == 0) stop("No variants overlap the supplied ranges")
    set_max <- max(vapply(ol, length, integer(1)))

    if (missing(title)) title <- paste0(
        scales::comma(length(ol_ids), 1), " Overlap; ",
        scales::comma(length(setdiff(ids, ol_ids)), 1),
        " No Overlap (", mcol, ")"
    )

    df_list <- list(ids = ol_ids, lapply(ol, \(x) ol_ids %in% x))
    df <- as.data.frame(df_list)

    if (!is.null(fill)) {
        fill <- match.arg(fill, colnames(mcols(gr)))
        extra_tbl <- data.frame(ids = mcols(gr)[[mcol]])
        extra_tbl[[fill]] <- mcols(gr)[[fill]]
        extra_tbl <- unique.data.frame(extra_tbl)
        i <- match(df$ids, extra_tbl$ids)
        df[[fill]] <- extra_tbl[[fill]][i]

    }

    sets <- c("Deletion", "Insertion", "SNV")
    stopifnot(is(fill_scale, "ScaleDiscrete"))
    p <- SimpleUpset::simpleUpSet(
        df, sets,
        set_layers = SimpleUpset::default_set_layers(
            fill = fill, hjust = hj_sets, expand = c(expand_sets, 0),
            label_size = label_size, fill_scale, guides(fill = guide_none())
        ),
        intersect_layers = SimpleUpset::default_intersect_layers(
            fill = fill, vjust = vj_intersect, expand = c(0, expand_intersect),
            label_size = label_size, fill_scale
        )
    ) + plot_annotation(title = title)
    p


}


