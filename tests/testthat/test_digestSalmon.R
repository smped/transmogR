f <- system.file("extdata/salmon_test", package = "transmogR")

test_that("assayFromQuants returns correct values",{

    ids <- c("t1", "t2")
    quants <- list(
        A = data.frame(Name = c("t1", "t2"), od = rnorm(2)),
        B = data.frame(Name = "t1", od = rnorm(1))
    )
    mat <- .assayFromQuants(quants, "od", ids, 0)
    expect_true(is(mat, "matrix"))
    expect_true(is.double(mat))
    expect_equal(rownames(mat), c("t1", "t2"))
    expect_equal(colnames(mat), c("A", "B"))
    expect_equal(mat[4], 0)

})

test_that("salmon digestion is smooth", {
    se <- suppressMessages(
        digestSalmon(f, extra_assays = c("TPM", "effectiveLength", "length"))
    )
    expect_true(is(se, "SummarizedExperiment"))
    expect_equal(dim(se), c(2L, 1L))
    expect_equal(metadata(se)$resampleType, "gibbs")
    expect_equal(
        c("counts", "scaledCounts", "TPM", "effectiveLength", "length"),
        assayNames(se)
    )
    expect_true(!any(is.na(rowData(se)$overdispersion)))

    se <- suppressMessages(digestSalmon(f))
    expect_equal(c("counts", "scaledCounts"), assayNames(se))

    se <- suppressMessages(digestSalmon(f, max_boot = 0))
    expect_equal("counts", assayNames(se))
    expect_equal("length", colnames(rowData(se)))

})

test_that("errors on incorrect directory", {
    expect_error(digestSalmon(f, aux_dir = "aux", verbose = FALSE))
    expect_error(digestSalmon(dirname(f), verbose = FALSE), "Missing json.+")
    expect_error(digestSalmon(file.path(f, "not_here")), "Unable.+")
})

test_that("errors on length_as_assay", {
    expect_warning(digestSalmon(f, length_as_assay = TRUE))
})

test_that("errors on single bootstrap", {
    expect_error(digestSalmon(f, max_boot = 1))
})

test_that("C parsing is correct",{
    f_nm <- file.path(f, "aux_info", "bootstrap", "names.tsv.gz")
    nm <- .Call("parse_trans_names", f_nm)
    nm_true <- c("ENST00000000233.10", "ENST00000000412.8")
    expect_equal(nm, nm_true)

    f_boot <- file.path(f, "aux_info", "bootstrap", "bootstraps.gz")
    od <- .C(
        "calc_boot_row_vals",
        filename = f_boot, n_trans = length(nm_true), n_boot = 1L,
        result = numeric(length(nm_true))
    )$result
    expect_equal(c(0, 0), od)

})
