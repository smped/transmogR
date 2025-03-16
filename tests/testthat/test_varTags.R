exons <- subset(gtf, type == "exon")
exon_by_trans <- splitAsList(exons, exons$transcript_id)[31:40]

test_that("tags are correct when produced", {

    tags <- varTags(exon_by_trans, var, tag = "mod")
    expect_true(length(tags) == length(exon_by_trans))
    expect_true(all(tags[c(2, 3, 5, 7:9)] == ""))
    expect_true(all(tags[c(1, 6, 10)] == "_mod_s"))
    expect_true(all(tags[c(4)] == "_mod_si"))

    empty <- varTags(exon_by_trans, var, var_tags = FALSE)
    expect_true(all(empty == ""))

    tag_only <- varTags(exon_by_trans, var, tag = "mod", var_tags = FALSE)
    expect_true(all(tag_only[c(1,4,6,10)] == "mod"))
    expect_true(all(tag_only[-c(1,4,6,10)] == ""))

})
