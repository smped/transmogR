library(GenomicFeatures)
library(GenomicRanges)
test_that("mogrifyTranscriptome works as expected",{
    if (requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) {
        hg19 <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
        id <- "ENST00000684155.1"
        ex <- GRanges(
            seqnames = "chrX",
            ranges = IRanges(
                start = c(49121123L, 49114753L, 49114121L, 49114002L),
                width = c(50L, 232L, 36L, 21L)
            ),
            strand  = "-", transcript_id = id, seqinfo = seqinfo(hg19)
        )
        grl <- GRangesList(ex)
        names(grl) <- id
        orig <- extractTranscriptSeqs(hg19, grl)
        var <- GRanges(
            c("chrX:49121172", "chrX:49114002-49114008"), seqinfo = seqinfo(hg19)
        )
        var$ALT <- c("G", "A") # SNP, Deletion at 1st & last positions
        new <- suppressMessages(
            mogrifyTranscriptome(
                hg19, var = var, exons = ex, tag = "test", var_tags = TRUE
            )
        )
        expect_true(names(new) == paste0(id, "_test_sd"))
        expect_true(as.character(new[[1]][1]) == "C")
        expect_true(length(new[[1]]) == 333)
        expect_true(as.character(new[[1]][seq(330, 333)]) == "GGAT")

        ## Check that omit_ranges filters unwanted regions
        empty <- suppressMessages(
            mogrifyTranscriptome(
                hg19, var = var, exons = ex, tag = "test", var_tags = TRUE,
                omit_ranges = range(ex)
            )
        )
        expect_true(length(empty) == 0)

        ## Check for an unaltered sequence if no variants
        new <- suppressMessages(
            mogrifyTranscriptome(hg19, var = var[NULL], exons = ex, tag = "test", var_tags = TRUE)
        )
        expect_equal(new, orig)

        ## Check for incorrect exons
        expect_error(
            suppressMessages(mogrifyTranscriptome(hg19,var = var[NULL], exons = unstrand(ex))),
            "Unstranded.+"
        )
    }
})


