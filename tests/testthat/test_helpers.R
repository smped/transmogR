test_that(".checkAlts behaves correctly", {
    gr <- GenomicRanges::GRanges("chr1:1")
    gr$ALT <- ""
    expect_error(.checkAlts(gr, "ALT"), "Please set Deletions.+")
    gr$ALT <- "x"
    expect_error(.checkAlts(gr, "ALT"), "Non-IUPAC.+")

    gr$ALT <- NA
    expect_error(
        suppressMessages(.checkAlts(gr, "ALT")), "All NA values.+"
    )

    gr$ALT <- "A"
    expect_equal(gr, .checkAlts(gr, "ALT"))
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

test_that(".makeIntersectionArgs behaves correctly",{
    new <- .makeIntersectionArgs(list(bar_number_threshold = 1, fill = "blue"))
    expect_true(is(new, "list"))
    types <- c(
        mapping = "call", counts = "logical", bar_number_threshold = "numeric",
        text_colors = "call", text = "call", text_mapping = "call", mode = "character",
        position = "call", width = "numeric", fill = "character"
    )
    expect_equal(vapply(new, \(x) is(x)[[1]], character(1)), types)
    expect_true(new$bar_number_threshold == 1)
    expect_true(new$fill == "blue")
    expect_error(.makeIntersectionArgs(NULL))
})
