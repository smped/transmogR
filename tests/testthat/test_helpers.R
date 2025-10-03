test_that(".checkAlts behaves correctly", {
    gr <- GenomicRanges::GRanges("chr1:1")
    gr$ALT <- ""
    expect_error(.checkAlts(gr, "ALT", ol_vars = "fail"), ".+Please set Deletions.+IUPAC.+")

    gr$ALT <- NA
    expect_error(.checkAlts(gr, "ALT", ol_vars = "fail"), ".+NA values.+")

    gr$ALT <- "A"
    expect_equal(gr, .checkAlts(gr, "ALT", ol_vars = "fail"))

    gr <- c(gr, gr)
    expect_error(.checkAlts(gr, "ALT", ol_vars = "fail"), ".+Duplicate.+")

})

test_that(".parseVariants behaves correctly",{
    library(VariantAnnotation)
    f <- VcfFile(
        system.file("extdata/1000GP_subset.vcf.gz", package = "transmogR")
    )
    gr <- .parseVariants(f, "ALT", GRanges("chr1:1-500000"))
    expect_equal(length(gr), 14L)
    expect_true(all(colnames(mcols(gr)) == c("REF", "ALT")))
    expect_true(is(gr$REF, "character"))
    expect_true(is(gr$ALT, "character"))

})


test_that(".checkOverlapVars works as expected", {
    test_var <- data.frame(
        seqnames = "chr10",
        start = c(2569741L, 2569750L, 12115100L, 12115100L, 12115104L),
        end = c(2569750L, 2569750L, 12115100L, 12115104L, 12115104L),
        strand = "*",
        REF = c("TTGTGTACTC", "C", "G", "GAACA", "A"),
        ALT = c("T", "CAAA", "C", "G", "ATTT")
    ) |> makeGRangesFromDataFrame(keep.extra.columns = TRUE)

    expect_error(
        .checkOverlapVars(test_var, "fail"), "5 pairs of overlapping loci.+"
    )
    expect_warning(
        gr <- .checkOverlapVars(test_var, "none"), ".+All overlapping.+"
    )
    expect_true(length(gr) == 0)

    expect_warning(
        gr <- .checkOverlapVars(test_var, "first"), ".+The first overlapping.+"
    )
    expect_equal(test_var[c(1, 3)], gr)

    expect_warning(
        gr <- .checkOverlapVars(test_var, "last"), ".+The last overlapping.+"
    )
    expect_equal(test_var[c(2, 5)], gr)

    expect_warning(
        gr <- .checkOverlapVars(test_var, "longest"), ".+The longest overlapping.+"
    )
    expect_equal(test_var[c(1, 4)], gr)

    expect_warning(
        gr <- .checkOverlapVars(test_var, "shortest"), ".+The shortest overlapping.+"
    )
    expect_equal(test_var[c(2, 3)], gr)


})

