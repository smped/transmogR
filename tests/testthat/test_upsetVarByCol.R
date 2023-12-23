test_that("upsetVarByCol behaves correctly", {
    grl <- splitAsList(gtf, gtf$type)
    p <- upsetVarByCol(grl$exon, var, intersection_args = list(fill = "blue"))
    expect_true(is(p, "patchwork"))
    expect_true(is(p[[1]], "spacer"))
    expect_true(is(p[[2]], "gg"))
    expect_true(p[[2]]$geom[[1]]$aes_params$fill == "blue")
    expect_true(is(p[[2]]$geom[[1]]$geom, "GeomBar"))
    expect_true(p[[2]]$labels$title == "37 Overlap; 75 No Overlap (transcript_id)")
    expect_true(is(p[[3]]$layers[[3]]$geom, "GeomText"))

    ## Turn off counts
    p <- upsetVarByCol(grl$exon, var, set_counts = FALSE)
    expect_true(length(p[[3]]$layers) == 2)

    ## Errors
    expect_error(upsetVarByCol(grl))
    expect_error(upsetVarByCol(grl$exon, ""))
    expect_error(upsetVarByCol(grl$exon, var, ""))

})
