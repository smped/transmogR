
test_that("varTypes calls variants correctly",{
    type <- varTypes(var)
    tb <- table(type)
    tb_true <- structure(
        c(Deletion = 6L, Insertion = 1L, SNV = 93L), dim = 3L,
        dimnames = list(type = c("Deletion", "Insertion", "SNV")),
        class = "table"
    )
    expect_equal(tb, tb_true)
})

test_that("varTypes errors as expected", {
    type <- varTypes(var)
    width(var)[type == "Insertion"] <- 2
    expect_error(varTypes(var), "Some variants.+")
    var$ALT <- 1
    expect_error(varTypes(var), ".+is not TRUE")
})
