f <- system.file("extdata/salmon_test", package = "transmogR")

test_that("assayFromQuants returns correct values",{

    ids <- c("t1", "t2")
    quants <- list(
        A = data.frame(Name = c("t1", "t2"), od = rnorm(2)),
        B = data.frame(Name = "t1", od = rnorm(1))
    )
    mat <- .assayFromQuants(quants, "od", ids, 0)
    expect_true(is(mat, "matrix"))
    expect_true(is.double(mat))
    expect_equal(rownames(mat), c("t1", "t2"))
    expect_equal(colnames(mat), c("A", "B"))
    expect_equal(mat[4], 0)

})

test_that("salmon digestion is smooth", {
 se <- suppressMessages(digestSalmon(f, length_as_assay = TRUE))
 expect_true(is(se, "SummarizedExperiment"))
 expect_equal(dim(se), c(2L, 1L))
 expect_equal(metadata(se)$resampleType, "gibbs")
 expect_equal(
     c("counts", "scaledCounts", "TPM", "effectiveLength", "length"),
     SummarizedExperiment::assayNames(se)
 )

 se <- suppressMessages(digestSalmon(f, extra_assays = NULL))
 expect_equal(
     c("counts", "scaledCounts"), SummarizedExperiment::assayNames(se)
 )

})

test_that("errors on incorrect directory", {
    expect_error(digestSalmon(f, aux_dir = "aux", verbose = FALSE))
    expect_error(digestSalmon(dirname(f), verbose = FALSE), "Missing json.+")
    expect_error(digestSalmon(file.path(f, "not_here")), "Unable.+")
})
