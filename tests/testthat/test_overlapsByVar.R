test_that("overlapsByVar behaves correctly", {

    grl <- splitAsList(gtf, gtf$type)
    ## Expected outputs
    mat <- overlapsByVar(grl, var)
    expect_true(is(mat, "matrix"))
    expect_equal(rownames(mat), names(grl))
    expect_equal(colnames(mat), c("Deletion", "Insertion", "SNV", "Total"))
    expect_true(all(mat[,"Total"] == c(98, 98, 38)))

    ## And the errors
    expect_error(overlapsByVar(GRanges()))
    expect_error(overlapsByVar(grl, ""))
    expect_error(overlapsByVar(grl, var, ""))
})
