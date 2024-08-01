

test_that("assayFromQuants returns correct values",{

    quants <- list(
        A = data.frame(Name = c("t1", "t2"), od = rnorm(2)),
        B = data.frame(Name = "t1", od = rnorm(1))
    )
    mat <- .assayFromQuants(quants, "od", 0)
    expect_true(is(mat, "matrix"))
    expect_true(is.double(mat))
    expect_equal(rownames(mat), c("t1", "t2"))
    expect_equal(colnames(mat), c("A", "B"))
    expect_equal(mat[4], 0)

})

test_that("salmon digestion is smooth", {
 se <- suppressMessages(digestSalmon(system.file("extdata", package = "transmogR")))
 expect_true(is(se, "SummarizedExperiment"))
 expect_equal(dim(se), c(2L, 1L))
 expect_equal(metadata(se)$resampleType, "gibbs")

})
