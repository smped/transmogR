library(GenomicRanges)
test_that("transmogrify works as expected",{
    if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
        hg38 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
        id <- "ENST00000684155.1"
        ex <- GRanges(
            seqnames = "chrX",
            ranges = IRanges(
                start = c(49264661L, 49258296L, 49257664L, 49257545L),
                width = c(50L, 232L, 36L, 21L)
            ),
            strand  = "-", transcript_id = id, seqinfo = seqinfo(hg38)
        )
        grl <- GRangesList(ex)
        names(grl) <- id
        orig <- GenomicFeatures::extractTranscriptSeqs(hg38, grl)
        var <- GRanges(
            c("chrX:49264710", "chrX:49257545-49257551"), seqinfo = seqinfo(hg38)
        )
        var$ALT <- c("G", "A") # SNP, Deletion at 1st & last positions
        new <- suppressMessages(
            transmogrify(
                hg38, var = var, exons = ex, tag = "test", var_tags = TRUE
            )
        )
        expect_true(names(new) == paste0(id, "_test_sd"))
        expect_true(as.character(new[[1]][1]) == "C")
        expect_true(length(new[[1]]) == 333)
        expect_true(as.character(new[[1]][seq(330, 333)]) == "GGAT")

        ## Check that omit_ranges filters unwanted regions
        empty <- suppressMessages(
            transmogrify(
                hg38, var = var, exons = ex, tag = "test", var_tags = TRUE,
                omit_ranges = range(ex)
            )
        )
        expect_true(length(empty) == 0)

        ## Check for an unaltered sequence if no variants
        new <- suppressMessages(
            transmogrify(hg38, var = var[NULL], exons = ex, tag = "test", var_tags = TRUE)
        )
        expect_equal(new, orig)

        ## Check for incorrect exons
        expect_error(
            suppressMessages(transmogrify(hg38,var = var[NULL], exons = unstrand(ex))),
            "Unstranded.+"
        )

        ## Now using a VCF
        exons <- rtracklayer::import.gff(
            system.file("extdata/gencode.v44.subset.gtf.gz", package = "transmogR"),
            feature.type = "exon"
        )
        vcf <- VariantAnnotation::VcfFile(
            system.file("extdata/1000GP_subset.vcf.gz", package = "transmogR")
        )
        new1 <- transmogrify(
            hg38, vcf, exons, tag = "test", verbose = FALSE,
            which = GRanges("chr1:100000-400000")
        )
        expect_true(sum(grepl("test", names(new1))) == 1)

        seq <- as(getSeq(hg38, names = "chr1"), "DNAStringSet")
        names(seq) <- "chr1"
        new2 <- transmogrify(
            seq, vcf, exons, tag = "test", verbose = FALSE,
            which = GRanges("chr1:100000-400000")
        )
        expect_equal(new1, new2)

    }
})


