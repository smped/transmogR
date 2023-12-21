library(VariantAnnotation)
f <- system.file("extdata/1000GP_subset.vcf.gz", package = "transmogR")
vcf <- readVcf(f)
gr <- rowRanges(vcf)

test_that("varTypes calls variants correctly",{
    type <- varTypes(gr)
    tb <- table(type)
    tb_true <- structure(
        c(Deletion = 6L, Insertion = 1L, SNV = 93L), dim = 3L,
        dimnames = list(type = c("Deletion", "Insertion", "SNV")),
        class = "table"
    )
    expect_equal(tb, tb_true)
})

test_that("varTypes errors as expected", {
    type <- varTypes(gr)
    width(gr)[type == "Insertion"] <- 2
    expect_error(varTypes(gr), "Some variants.+")
    gr$ALT <- 1
    expect_error(varTypes(gr), ".+is not TRUE")
})
