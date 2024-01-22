test_that("overlapsByVar behaves correctly", {

    grl <- splitAsList(gtf, gtf$type)
    ## Expected outputs
    df <- overlapsByVar(grl, var)
    expect_true(is(df, "data.frame"))
    expect_equal(rownames(df), names(grl))
    expect_equal(colnames(df), c("Deletion", "Insertion", "SNV", "Total"))
    expect_true(all(df$Total == c(98, 98, 38)))

    ## And the errors
    expect_error(overlapsByVar(GRanges()))
    expect_error(overlapsByVar(grl, ""))
    expect_error(overlapsByVar(grl, var, ""))
})
