library(VariantAnnotation)
f <- system.file("extdata/1000GP_subset.vcf.gz", package = "transmogR")
vcf <- readVcf(f)
gr <- rowRanges(vcf)

test_that("calvInDel calls variants correctly",{
    type <- calvInDel(gr)
    tb <- table(type)
    tb_true <- structure(
        c(Deletion = 6L, Insertion = 1L, SNV = 93L), dim = 3L,
        dimnames = list(type = c("Deletion", "Insertion", "SNV")),
        class = "table"
    )
    expect_equal(tb, tb_true)
})

test_that("calvInDel errors as expected", {
    type <- calvInDel(gr)
    width(gr)[type == "Insertion"] <- 2
    expect_error(calvInDel(gr), "Some variants.+")
    gr$ALT <- 1
    expect_error(calvInDel(gr), ".+is not TRUE")
})
