test_that(".checkAlts functions correctly", {
    gr <- GenomicRanges::GRanges("chr1:1")
    gr$ALT <- ""
    expect_error(.checkAlts(gr, "ALT"), "Non-IUPAC.+")

    gr$ALT <- NA
    expect_message(.checkAlts(gr, "ALT"), "NA values.+")

    gr$ALT <- "A"
    expect_equal(gr, .checkAlts(gr, "ALT"))
})
