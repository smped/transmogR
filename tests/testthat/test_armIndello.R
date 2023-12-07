library(GenomicRanges)
test_that("armIndello replaces InDels correctly for a DNAString", {
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
    result <- armIndello(seq, indels, exons, alt_col = "ALT")

    ## Check if the replacements are correct
    expect_equal(as.character(result), "AAACGATCCCG")

    ## Test the reverse strand
    strand(exons) <- "-"
    result <- armIndello(seq, indels, exons, alt_col = "ALT")
    expect_equal(as.character(result), "ATGTCGATTTG")

    ## Check the return for empty indels
    result <- armIndello(seq, indels[NULL], exons, alt_col = "ALT")
    expect_equal(as.character(result), as.character(seq))

})

test_that("armIndello replaces Indels for a DNAStringSet", {

    ## Create a sample DNA sequence
    seq <- DNAStringSet(c(seq1 = "ATCGATCGATCG"))
    ## Create a sample GenomicRanges object with indels
    indels <- GRanges(seqnames = "seq1",
                      ranges = IRanges(start = c(2, 8), end = c(2, 10)),
                      ALT = c("AA", "C"))
    result <- suppressMessages(armIndello(seq, indels))
    expect_true(length(result) == 1)
    expect_true(is(result, "DNAStringSet"))
    expect_equal(as.character(result[["seq1"]]), "AAACGATCCCG")

    ## Check output when no indels
    expect_message(armIndello(seq, indels[NULL]), "No InDels.+")

})

test_that("armIndello replaces InDels for BSgenome", {
    run_test <- requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
    if (run_test) {
        bs <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
        indels <- GRanges("chrM:1-3", seqinfo = seqinfo(bs))
        indels$ALT <- "G"
        result <- suppressMessages(
            armIndello(bs, indels, names = "chrM", alt_col = "ALT")
        )
        expect_true(length(result[["chrM"]]) == 16569)
        expect_equal(as.character(result[["chrM"]][1:3]), "GCA")

    }
})
