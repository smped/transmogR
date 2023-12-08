library(GenomicRanges)
test_that("mogrifyGenome works correctly", {
    if (requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) {
        hg19 <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
        var <- GRanges(c("chrM:1", "chrM:2"), seqinfo = seqinfo(hg19))
        var$ALT <- c("A", "NN") # Normally a G
        new <- suppressMessages(
            mogrifyGenome(hg19, var, alt_col = "ALT", names = "chrM", tag = "test", var_tags = TRUE)
        )
        expect_true(as.character(subseq(new[[1]], 1, 3)) == "ANN")
        expect_true(names(new)[[1]] == "chrM_test_si")

    }
    # Test for an empty set of variants
    seq <- DNAStringSet(c(seq1 = "AAT"))
    var <- GRanges("seq1:1", seqinfo = seqinfo(seq))
    var$ALT <- "A"
    unch <- mogrifyGenome(seq, var[NULL], tag = "test")
    expect_true(length(unch[[1]]) == 3)
    expect_true(names(unch)[[1]] == "seq1")

})
