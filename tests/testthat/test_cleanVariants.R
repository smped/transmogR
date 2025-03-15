test_that("Parsing works", {
  var <- cleanVariants(VcfFile(vcf))
  expect_true(is(var, "GRanges"))
  expect_true(all(vapply(mcols(var), \(x) is.character(x), logical(1))))
})
