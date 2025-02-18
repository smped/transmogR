#' Parse the output from salmon
#'
#' Parse transcript counts and additional data from salmon
#'
#' @details
#' This function is based heavily on [edgeR::catchSalmon()] however, there are
#' some important differences:
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
#' This greater flexibility also requires more stringent checking and, as such,
#' for smaller datasets, digestSalmon may be slower that the edgeR function.
#'
#' The SummarizedExperiment object returned may also contain multiple assays,
#' as described elsewhere on this page
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
#' @param max_boot The maximum number of bootstraps to use. Setting this to
#' zero will ignore all bootstraps and the scaledCounts assay will not be
#' included in the returned object
#' @param n_threads The number of threads to use for all parallel processes
#' @param c_threads Passed to C for optional parallelised computation of
#' bootstrap overdispersion estimates. Only used if OpenMP is available on
#' the underlying computational infrastructure.
#' @param r_threads Number of threads used for any internal calls to
#' [parallel::mclapply()]. Given the default memory management practices in R,
#' setting this value too high may cause issues with RAM requirements
#' @param ... Not used
#'
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame metadata<-
#' @importFrom matrixStats rowVars
#'
#' @export
digestSalmon <- function(
        paths, max_sets = 2L, aux_dir = "aux_info", name_fun = basename,
        verbose = TRUE, extra_assays = NULL, max_boot = Inf, n_threads = 1,
        ..., r_threads = n_threads, c_threads = n_threads
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
    n_boot <- as.integer(min(n_boot, max_boot))
    if (n_boot == 1) stop("The number of bootstraps cannot be equal to 1")
    r_threads <- as.integer(r_threads)
    c_threads <- as.integer(c_threads)

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

    ## Column Types to ensure only the required columns are parsed
    col_types <- list(
        Name = "c", Length = "d", effectiveLength = "-", TPM = "-",
        NumReads = "d"
    )
    col_types[names(col_types) %in% extra_assays] <- "d"
    col_types <- paste(unlist(col_types), collapse = "")

    ## Import quants
    options(readr.show_progress = FALSE)
    if (verbose) message("Parsing quants...")
    quants <- lapply(quant_files, vroom::vroom, col_types = col_types)
    if (verbose) message("done")

    ## Transcript Lengths
    if (verbose) message("Checking transcript lengths...")
    ids <- sort(unique(unlist(lapply(quants, \(x) x$Name))))
    trans_len <- .assayFromQuants(quants, "Length", ids, NA_integer_, r_threads)
    if (!("length" %in% extra_assays)) {
        if (length(quants) > 1 & any(rowVars(trans_len, na.rm = TRUE) > 0)) {
            msg <- paste(
                "Some transcripts have differing lengths between samples.",
                "Please set extra_assays = 'length'"
            )
            stop(msg)
        }
        # Delete the object for downstream if not required
        trans_len <- NULL
    }
    if (verbose) message("done")

    ## Setup the rowData & assays
    if (verbose) message("Obtaining assays...")
    counts <- .assayFromQuants(quants, "NumReads", ids, 0, r_threads)
    assays <- list(counts = counts)
    rowDF <- DataFrame(row.names = ids)
    md <- list(resampleType = boot_types)
    if (n_boot > 0) {
        if (verbose) message("Estimating overdispersions...")
        final_od <- .overdispFromBoots(
            paths, n_boot, .ids = ids, n_threads = c_threads
        )
        if (verbose) message("done")
        assays$scaledCounts <- counts / final_od
        rowDF$overdispersion <- final_od
        md$n_boot <- n_boot
    }
    ## Extra Assays
    if ("TPM" %in% extra_assays)
        assays$TPM <- .assayFromQuants(quants, "TPM", ids, 0, r_threads)
    if ("effectiveLength" %in% extra_assays)
        assays$effectiveLength <- .assayFromQuants(
            quants, "EffectiveLength", ids, NA_real_, r_threads
        )
    assays$length <- trans_len
    if (verbose) message("done")

    ## Handle a single sample case where R defaults to vectors
    if (length(paths) == 1) assays <- lapply(assays, as.matrix)
    if (!("length" %in% extra_assays)) {
        i <- which.max(n_trans)
        rowDF$length <- setNames(quants[[i]]$Length, quants[[i]]$Name)[ids]
    }

    colDF <- DataFrame(totals = colSums(assays$counts), n_trans = n_trans)
    se <- SummarizedExperiment(assays = assays, rowData = rowDF, colData = colDF)
    metadata(se) <- md
    colnames(se) <- paths

    if (is(name_fun, "function")) colnames(se) <- name_fun(colnames(se))
    se

}


#' @useDynLib transmogR, .registration = TRUE
#' @importFrom matrixStats rowMeans2 rowSums2
#' @importFrom parallel mclapply
#' @importFrom stats setNames median qf
#' @keywords internal
.overdispFromBoots <- function(paths, n_boot, .ids, n_threads) {

    suf <- file.path("aux_info", "bootstrap", "bootstraps.gz")
    boot_files <- vapply(paths, file.path, character(1), suf)
    if (!all(file.exists(boot_files))) stop("Missing bootstrap files")
    id_files <- gsub("bootstraps.gz$", "names.tsv.gz", boot_files)
    if (!all(file.exists(id_files))) stop("Missing names.tsv.gz files")

    ## Try a more computationally efficient approach
    sums_ti <- lapply(
        seq_along(paths),
        \(i){
            trans_ids <- .Call("parse_trans_names", id_files[[i]])
            n <- length(trans_ids)
            f <- boot_files[[i]]
            sum_ti <- .Call("calc_boot_row_vals", f, n, n_boot, n_threads)
            names(sum_ti) <- trans_ids
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

.assayFromQuants <- function(x, var, .ids, fill = NA_real_, r_threads) {

    mat <- do.call(
        "cbind",
        mclapply(
            x, \(x) setNames(x[[var]], x[["Name"]])[.ids], mc.cores = r_threads
        )
    )
    mat[is.na(mat)] <- fill
    mat

}


