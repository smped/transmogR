test_that("upsetVarByCol behaves correctly", {
    grl <- splitAsList(gtf, gtf$type)
    p <- upsetVarByCol(
        grl$exon, var, fill = "gene_type", fill_scale = scale_fill_brewer(palette = "Set1")
    )
    expect_true(is(p, "patchwork"))
    expect_true(is_waiver(p[[1]]@data))
    expect_equal(
        p@meta$patches$annotation$title, "37 Overlap; 75 No Overlap (transcript_id)"
    )
    expect_equal(
        as_label(p[[2]]@layers$geom_bar$mapping$fill), "gene_type"
    )
    expect_equal(
        as_label(p[[3]]@layers$geom_bar$mapping$fill), "gene_type"
    )
    expect_true(is(p[[2]]@scales$scales[[3]], "ScaleDiscrete"))
    expect_equal(p[[2]]@scales$scales[[3]]$aesthetics, "fill")
    expect_equal(
        p[[2]]@scales$scales[[3]]$palette(3), c("#E41A1C", "#377EB8", "#4DAF4A")
    )
    expect_true(is(p[[3]]@guides$guides$fill, "GuideNone"))

    ## Errors
    expect_error(upsetVarByCol(grl))
    expect_error(upsetVarByCol(grl$exon, ""))
    expect_error(upsetVarByCol(grl$exon, var, ""))
    expect_error(upsetVarByCol(grl$exon, var[NULL]), "No variants overlap.+")

})
