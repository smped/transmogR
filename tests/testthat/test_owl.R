test_that("owl substitutes SNPs correctly", {
    # Create a sample XStringSet
    seq <- DNAStringSet(c("ATCG", "GATC", "CGTA"))
    names(seq) <- paste0("seq", 1:3)

    # Create a sample GenomicRanges object with SNPs
    snps <- GRanges(
        seqnames = c("seq1", "seq2", "seq3"),
        ranges = IRanges(start = c(2, 1, 3), end = c(2, 1, 3)),
        ALT = c("C", "T", "A")
    )

    # Call the owl function
    result <- owl(seq, snps, alt_col = "ALT")

    # Check if the substitutions are correct
    expect_equal(
        as.character(result), c(seq1 = "ACCG", seq2 = "TATC", seq3 = "CGAA")
    )
})

test_that("owl errors when expected", {
    # Create a sample XStringSet
    seq <- DNAStringSet(c("ATCG", "GATC", "CGTA"))

    # Create a sample GenomicRanges object with SNPs
    snps <- GRanges(
        seqnames = c("seq1", "seq2", "seq3"),
        ranges = IRanges(start = c(2, 1, 3), end = c(2, 1, 3)),
        ALT = c("C", "T", "A")
    )

    expect_error(owl(seq, snps, alt_col = "ALT"), "all\\(seqlevels\\(.+TRUE$")
    expect_error(owl("", snps, alt_col = "ALT"))
    expect_error(owl(seq, "", alt_col = "ALT"))

})
