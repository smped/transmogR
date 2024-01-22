test_that("parY works & errors correctly", {

    # The correct length for Y
    sq <- GenomeInfoDb::Seqinfo(
        seqnames = "chrY", seqlengths = 59373566, genome = "hg19_only_chrY"
    )
    parY <- parY(sq)
    expect_equal(start(parY), c(10001, 59034050))

    sq <- GenomeInfoDb::Seqinfo(
        seqnames = "chrY", seqlengths = 5937366, genome = "incorrect_chrY"
    )
    expect_error(parY(sq), "Invalid.+")

    hg19 <- GRanges(c("Y:10001-2649520", "Y:59034050-59363566"))
    hg19$PAR <- paste0("PAR", 1:2)
    genome(hg19) <- "hg19"
    expect_equal(parY("hg19"), hg19)

})
