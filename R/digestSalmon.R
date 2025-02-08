#' Parse the output from salmon
#'
#' Parse transcript counts and additional data from salmon
#'
#' @details
#' This function is based heavily on [edgeR::catchSalmon()] with some important
#' exceptions:
#'
#' 1. A SummarizedExperiment object is returned
#' 2. Differing numbers of transcripts are allowed between samples
#'
#' The second point is intended for the scenario where some samples may have
#' been aligned to a full reference, with remaining samples aligned to a
#' partially masked reference (e.g. chrY). This will lead to differing numbers
#' of transcripts within each salmon index, however, common estimates of
#' overdispersions are required for scaling transcript-level counts. By default,
#' the function will error if >2 different sets of transcripts are detected,
#' however this can be modified using the max_sets argument.
#'
#' The SummarizedExperiment object returned will also contain multiple assays,
#' as described below
#'
#' @return A SummarizedExperiment object containing assays for counts and
#' scaledCounts.
#' The scaledCounts assay contains counts divided by overdispersions.
#' rowData in the returned object will also include transcript-lengths along
#' with the overdispersion estimates used to return the scaled counts.
#' TPM, effectiveLength and length can be returned as additional assays by
#' specifying one or more of these in the extra_assays argument
#'
#' @param paths Vector of file paths to directories containing salmon results
#' @param max_sets The maximum number of indexes permitted
#' @param aux_dir Subdirectory where bootstraps and meta_info.json are stored
#' @param name_fun Function applied to paths to provide colnames in the returned
#' object. Set to NULL or c() to disable.
#' @param verbose Print progress messages
#' @param extra_assays Can take values in  c("TPM", "effectiveLength", "length")
#' to optionally request TPM, effectiveLength or length as assays. Including
#' the length assay is intended for the use case of personalised transcriptomes
#' where transcript lengths may no longer be uniform across samples.
#' None will be returned by default
#' @param max_boot The maximum number of bootstraps to use
#' @param ... Not used
#'
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame metadata<-
#'
#' @export
digestSalmon <- function(
        paths, max_sets = 2L, aux_dir = "aux_info", name_fun = basename,
        verbose = TRUE, extra_assays = NULL, max_boot = Inf, ...
) {

    ## Initial file.path checks
    dir_exists <- vapply(paths, dir.exists, logical(1))
    if (!all(dir_exists)) {
        msg <- paste("Unable to find:", paths[!dir_exists], sep = "\n")
        stop(msg)
    }

    ## Handle the extra_assays & change in arguments
    if ("length_as_assay" %in% names(list(...))) {
        msg <- paste(
            "The argument 'length_as_assay' has been deprecated with v1.1.4.",
            "Please pass 'length' to the argument extra_assays.",
            sep = "\n"
        )
        warning(msg)
    }
    if (!is.null(extra_assays)) {
        valid_assays <- c("TPM", "effectiveLength", "length")
        extra_assays <- match.arg(extra_assays, valid_assays, several.ok = TRUE)
    }

    ## json checks
    if (verbose) message("Parsing json metadata...", appendLF = FALSE)
    ## Instead of getting aux_info from the cmd_info.json file, require this
    ## to be passed using the aux_dir argument
    aux_dir <- rep_len(aux_dir, length(paths))
    meta_json <- file.path(paths, aux_dir, "meta_info.json")
    json_exists <- file.exists(meta_json)
    if (!all(json_exists)) {
        msg <- paste("Missing json files:", meta_json[!json_exists], sep = "\n")
        stop(msg)
    }
    meta_info <- lapply(meta_json, jsonlite::fromJSON)
    ## Check bootstrap info
    boot_types <- unique(vapply(meta_info, \(x) x$samp_type, character(1)))
    if (length(boot_types) > 1) stop("Bootstraps must all use the same method")
    n_boot <- vapply(meta_info, \(x) max(x$num_bootstraps, 0L), integer(1))
    n_boot <- min(n_boot, max_boot)

    ## Check the transcriptomes
    n_trans <- vapply(
        meta_info,
        \(x) unlist(x[c("num_targets", "num_valid_targets")])[[1]],
        integer(1)
    )
    if (length(n_trans) != length(paths)) stop("Missing values in json files")
    n_sets <- length(unique(n_trans))
    if (n_sets > max_sets) stop(n_sets, " sets of annotations detected")
    if (verbose) message("done")

    ## quant checks
    quant_files <- vapply(paths, file.path, character(1), "quant.sf")
    quant_exists <- file.exists(quant_files)
    if (!all(quant_exists)) {
        msg <- paste("Missing quant files:", quant_files[!quant_exists], sep = "\n")
        stop(msg)
    }

    ## Import quants
    options(readr.show_progress = FALSE)
    if (verbose) message("Parsing quants...")
    quants <- lapply(quant_files, vroom::vroom, col_types = "cdddd")
    if (verbose) message("done")

    ## Transcript Lengths
    if (verbose) message("Checking transcript lengths...")
    lens <- unique(
        do.call("rbind", lapply(quants, \(x) x[c("Name", "Length")]))
    )
    ids <- lens[["Name"]]
    ## Handle transcripts which have multiple lengths, as may be the case for a
    ## set of personalised references
    if (any(duplicated(ids)) & !("length" %in% extra_assays)) {
        msg <- paste(
            "Some transcripts have differing lengths between samples.",
            "Please set extra_assays = 'length'"
        )
        stop(msg)
    }
    ## If passing here, ids will be unique or lengths will be an assay
    ## Setting as unique is still needed if passing to an assay
    ids <- unique(ids)
    if (verbose) message("done")

    # ## Setup the core assays
    if (verbose) message("Obtaining assays...")
    counts <- .assayFromQuants(quants, "NumReads", ids, 0)
    tpm <- eff_len <- trans_len <- NULL # Default to NULL
    if ("TPM" %in% extra_assays)
        tpm <- .assayFromQuants(quants, "TPM", ids, 0)
    if ("effectiveLength" %in% extra_assays)
        eff_len <- .assayFromQuants(quants, "EffectiveLength", ids, NA_real_)
    if ("length" %in% extra_assays)
        trans_len <- .assayFromQuants(quants, "Length", ids, NA_integer_)
    if (verbose) message("done")

    ## Now the bootstraps.
    final_od <- 1
    if (n_boot > 0) {
        if (verbose) message("Estimating overdispersions...")
        final_od <- .overdispFromBoots(paths, n_boot, n_trans, quants, .ids = ids)
        if (verbose) message("done")
    }

    ## All of the above can go into an SE.
    assays <- list(counts = counts, scaledCounts = counts / final_od)
    assays$TPM <- tpm
    assays$effectiveLength <- eff_len
    assays$length <- trans_len

    ## Handle a single sample case where R defaults to vectors
    if (length(paths) == 1) assays <- lapply(assays, as.matrix)
    rowDF <- DataFrame(overdispersion = final_od, row.names = ids)
    if (!("length" %in% extra_assays)) rowDF$length <- lens[["Length"]]

    colDF <- DataFrame(totals = colSums(assays$counts), n_trans = n_trans)
    se <- SummarizedExperiment(assays = assays, rowData = rowDF, colData = colDF)
    metadata(se) <- list(resampleType = boot_types)
    colnames(se) <- paths

    if (is(name_fun, "function")) colnames(se) <- name_fun(colnames(se))
    se

}

#' @importFrom matrixStats rowMeans2 rowSums2
#' @importFrom stats setNames median qf
#' @keywords internal
.overdispFromBoots <- function(paths, n_boot, n_trans, quants, .ids) {

    suf <- file.path("aux_info", "bootstrap", "bootstraps.gz")
    boot_files <- vapply(paths, file.path, character(1), suf)
    if (!all(file.exists(boot_files))) stop("Missing bootstrap files")

    ## Try a more computationally efficient approach
    sums_ti <- lapply(
        seq_along(paths),
        \(i){

            con <- gzcon(file(boot_files[[i]], open = "rb"))
            ## Enable different numbers of transcripts for different references
            boots <- readBin(con, what = 'double', n = n_trans[[i]] * n_boot)
            close(con)
            dim(boots) <- c(n_trans[[i]], n_boot)
            rownames(boots) <- quants[[i]]$Name
            lambda_ti <- rowMeans2(boots)
            sum_ti <- rowSums2((boots - lambda_ti)^2) / lambda_ti
            ## Return zero for undetectable transcripts (lambda_ti == 0)
            sum_ti[is.nan(sum_ti)] <- 0
            sum_ti[.ids]
        }
    )

    ## Following Baldoni et al, except d_t is not n(B - 1) where the transcript
    ## was not > 0 in all libraries
    boot_sums <- do.call("cbind", sums_ti)
    d_t <- rowSums2(boot_sums > 0, na.rm = TRUE) * (n_boot - 1)
    od_t <- rowSums2(boot_sums, na.rm = TRUE) / d_t
    ## Key values for the moderation of the overdispersion
    i <- d_t > 0
    d_0 <- 3
    d_med <- median(d_t[i])
    od_0 <- max(1, median(od_t[i]) / qf(0.5, d_med, d_0))
    ## Moderate the overdispersions
    od_mod <- pmax(1, (d_0 * od_0 + d_t * od_t) / (d_0 + d_t))
    od_mod[is.na(od_mod)] <- od_0
    od_mod

}

.assayFromQuants <- function(x, var, .ids, fill = NA_real_) {

    mat <- do.call(
        "cbind", lapply(x, \(x) setNames(x[[var]], x[["Name"]])[.ids])
    )
    mat[is.na(mat)] <- fill
    mat

}


