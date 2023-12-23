test_that("sjFromExons behaves correctly", {
    ex <- subset(gtf, type == "exon")
    gr <- sjFromExons(ex)
    expect_true(setdiff(colnames(mcols(gr)), colnames(mcols(ex))) == "site")
    expect_true(all(gr$site %in% c("donor", "acceptor")))
    expect_true(is(gr$transcript_id, "CharacterList"))
    expect_equal(length(gr), 196)

    gr <- sjFromExons(ex, chop = FALSE)
    expect_true(is(gr$transcript_id, "character"))
    expect_equal(length(gr), 730)

    ## Errors
    expect_error(sjFromExons(NULL))
    expect_error(sjFromExons(ex, rank_col = ""))
    expect_error(sjFromExons(ex, tx_col = ""))
    expect_error(sjFromExons(unstrand(ex)))
})
