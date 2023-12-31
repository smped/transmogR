library(GenomicRanges)
test_that("indelcator replaces InDels correctly for a DNAString", {
    ## Create a sample DNA sequence
    seq <- DNAString("ATCGATCGATCG")

    ## Create a sample GenomicRanges object with indels
    indels <- GRanges(seqnames = "seq1",
                      ranges = IRanges(start = c(2, 8), end = c(2, 10)),
                      ALT = c("AA", "C"))

    ## Create a sample GRanges object representing exons
    exons <- GRanges(
        seqnames = "seq1", ranges = IRanges(start = c(1, 7), end = c(6, 12))
    )

    ## Call the replaceVariants function
    result <- indelcator(seq, indels, exons, alt_col = "ALT")

    ## Check if the replacements are correct
    expect_equal(as.character(result), "AAACGATCCCG")

    ## Test the reverse strand
    strand(exons) <- "-"
    result <- indelcator(seq, indels, exons, alt_col = "ALT")
    expect_equal(as.character(result), "ATGTCGATTTG")

    ## Check the return for empty indels
    result <- indelcator(seq, indels[NULL], exons, alt_col = "ALT")
    expect_equal(as.character(result), as.character(seq))

})

test_that("indelcator replaces Indels for a DNAStringSet", {

    ## Create a sample DNA sequence
    seq <- DNAStringSet(c(seq1 = "ATCGATCGATCG"))
    ## Create a sample GenomicRanges object with indels
    indels <- GRanges(seqnames = "seq1",
                      ranges = IRanges(start = c(2, 8), end = c(2, 10)),
                      ALT = c("AA", "C"))
    result <- suppressMessages(indelcator(seq, indels))
    expect_true(length(result) == 1)
    expect_true(is(result, "DNAStringSet"))
    expect_equal(as.character(result[["seq1"]]), "AAACGATCCCG")

    ## Check output when no indels
    expect_message(indelcator(seq, indels[NULL]), "No InDels.+")

})

test_that("indelcator replaces InDels for BSgenome", {
    run_test <- requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)
    if (run_test) {
        bs <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
        indels <- GRanges("chrM:1-3", seqinfo = seqinfo(bs))
        indels$ALT <- "G"
        result <- suppressMessages(
            indelcator(bs, indels, names = "chrM", alt_col = "ALT")
        )
        expect_true(length(result[["chrM"]]) == 16567)
        expect_equal(as.character(result[["chrM"]][1:3]), "GCA")

    }
})
