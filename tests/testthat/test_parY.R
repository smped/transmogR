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

})
