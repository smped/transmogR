library(GenomicRanges)
library(VariantAnnotation)
test_that("genomogrify works correctly", {
    if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
        ## Using a GRanges with variants
        hg38 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
        var <- GRanges(c("chrM:1", "chrM:2"), seqinfo = seqinfo(hg38))
        var$ALT <- c("A", "NN") # Normally a G
        new <- genomogrify(
            hg38, var, alt_col = "ALT", names = "chrM", tag = "test",
            var_tags = TRUE, verbose = FALSE
        )
        expect_true(as.character(subseq(new[[1]], 1, 3)) == "ANN")
        expect_true(names(new)[[1]] == "chrM_test_si")

        ## Using a VCF
        vcf <- VcfFile(
            system.file("extdata/1000GP_subset.vcf.gz", package = "transmogR")
        )
        new <- genomogrify(
            hg38, vcf, names = "chr1", which = GRanges("chr1:839500-839550"),
            verbose = FALSE
        ) # Just one deletion
        expect_true(length(new[["chr1"]]) == 248956394)

        ## Check the masking
        new <- genomogrify(
            hg38, var, alt_col = "ALT", mask = GRanges("chrM:16561-16569"),
            names = "chrM", tag = "test", var_tags = TRUE, verbose = FALSE
        )
        expect_equal(as.character(new[[1]][16562:16570]), "NNNNNNNNN")

    }

    # Test for an empty set of variants
    seq <- DNAStringSet(c(seq1 = "AAT"))
    var <- GRanges("seq1:1", seqinfo = seqinfo(seq))
    var$ALT <- "A"
    unch <- genomogrify(seq, var[NULL], tag = "test")
    expect_true(length(unch[[1]]) == 3)
    expect_true(names(unch)[[1]] == "seq1")

    ## Test for an XStringSet
    test_seq <- as(paste(rep("N", 113970), collapse = ""), "DNAStringSet")
    names(test_seq) <- "chr1"
    new <- suppressMessages(
        genomogrify(test_seq, vcf, which = GRanges("chr1:1-113970"))
    )
    expect_true(as.character(new[["chr1"]][113969:113970]) == "TN")

})
