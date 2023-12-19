test_that(".checkAlts functions correctly", {
    gr <- GenomicRanges::GRanges("chr1:1")
    gr$ALT <- ""
    expect_error(.checkAlts(gr, "ALT"), "Non-IUPAC.+")

    gr$ALT <- NA
    expect_message(.checkAlts(gr, "ALT"), "NA values.+")

    gr$ALT <- "A"
    expect_equal(gr, .checkAlts(gr, "ALT"))
})

test_that(".pasreVariants functions correctly",{
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
