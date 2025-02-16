f <- system.file("extdata/salmon_test", package = "transmogR")
f_boot <- file.path(f, "aux_info", "bootstrap", "bootstraps.gz")
f_nm <- file.path(f, "aux_info", "bootstrap", "names.tsv.gz")

test_that("assayFromQuants returns correct values",{

    ids <- c("t1", "t2")
    quants <- list(
        A = data.frame(Name = c("t1", "t2"), od = rnorm(2)),
        B = data.frame(Name = "t1", od = rnorm(1))
    )
    mat <- .assayFromQuants(quants, "od", ids, 0, 1L)
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

    nm <- .Call("parse_trans_names", f_nm)
    nm_true <- c("ENST00000000233.10", "ENST00000000412.8")
    expect_equal(nm, nm_true)

    od <- .Call("calc_boot_row_vals", f_boot, length(nm_true), 10L, 1L)
    expect_true(all.equal(c(623.060848059843, 1813.5696136543), od))

})

test_that("C errors correctly when parsing transcript ids", {

    expect_error(.Call("parse_trans_names", ""), "Failed to open file.+")
    expect_error(
        .Call("parse_trans_names", NULL),
        "Filename must be a single character string"
    )
    expect_error(
        .Call("parse_trans_names", c(f_nm, f_nm)),
        "Filename must be a single character string"
    )
    expect_error(
        .Call("parse_trans_names", f_boot), "No tabs found in the file"
    )

    f_empty <- tempfile(fileext = "gz")
    f_con <- gzcon(file(f_empty, "wb"))
    cat(NULL, file = f_con)
    close(f_con)
    expect_error(.Call("parse_trans_names",  f_empty), "Empty file")
    unlink(f_empty)


})

test_that("C errors correctly when processing bootstraps", {
    expect_error(
        .Call("calc_boot_row_vals", f_boot, 1L, 1L, 1L), ".+transcripts.+"
    )
    expect_error(
        .Call("calc_boot_row_vals", f_boot, 2L, 1L, 1L), ".+bootstrap.+"
    )
    expect_error(
        .Call("calc_boot_row_vals", f_boot, 2L, 10L, 0L), ".+threads.+"
    )
    expect_error(
        .Call("calc_boot_row_vals", "", 2L, 10L, 1L), "Failed to open+"
    )
    expect_warning(
        .Call("calc_boot_row_vals", f_boot, 2L, 2L, 1L),
        "Additional data exists in file.+"
    )
})
