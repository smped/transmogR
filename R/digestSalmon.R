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
#' @return A SummarizedExperiment object containing assays for counts,
#' scaledCounts, TPM and effectiveLength.
#' The scaledCounts assay contains counts divided by overdispersions.
#' rowData in the returned object will also include transcript-lengths along
#' with the overdispersion estimates used to return the scaled counts.
#'
#' @param paths Vector of file paths to directories containing salmon results
#' @param max_sets The maximum number of indexes permitted
#' @param aux_dir Subdirectory where bootstraps and meta_info.json are stored
#' @param name_fun Function applied to paths to provide colnames in the returned
#' object. Set to NULL or c() to disable.
#' @param verbose Print progress messages
#'
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame metadata<-
#'
#' @export
digestSalmon <- function(
        paths, max_sets = 2L, aux_dir = "aux_info", name_fun = basename,
        verbose = TRUE
) {

    ## Initial file.path checks
    dir_exists <- vapply(paths, dir.exists, logical(1))
    if (!all(dir_exists)) {
        msg <- paste("Unable to find:", paths[!dir_exists], sep = "\n")
        stop(msg)
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
    n_trans <- vapply(
        meta_info,
        \(x) unlist(x[c("num_targets", "num_valid_targets")])[[1]],
        integer(1)
    )
    n_sets <- length(unique(n_trans))
    n_boot <- vapply(meta_info, \(x) max(x$num_bootstraps, 0L), integer(1))
    if (n_sets > max_sets) stop(n_sets, " sets of annotations detected")
    if (length(n_trans) != length(n_boot)) stop("Missing values in json files")
    types <- unique(vapply(meta_info, \(x) x$samp_type, character(1)))
    if (length(types) > 1) stop("Bootstraps must all use the same method")
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
    lens <- unique(do.call("rbind", quants)[c("Name", "Length")])
    ids <- lens[["Name"]]
    ## Setup the core assays
    if (verbose) message("Obtaining assays...", appendLF = FALSE)
    counts <- .assayFromQuants(quants, "NumReads", 0)[ids, ]
    tpm <- .assayFromQuants(quants, "TPM", 0)[ids, ]
    eff_len <- .assayFromQuants(quants, "EffectiveLength", NA_real_)[ids, ]
    if (verbose) message("done")

    ## Now the bootstraps.
    final_od <- 1
    if (any(n_boot > 0)) {
        if (any(n_boot == 0)) warning("Some samples were not bootstrapped")
        if (verbose) message("Estimating overdispersions...", appendLF = FALSE)
        final_od <- .overdispFromBoots(paths, n_boot, n_trans, quants)[ids]
        if (verbose) message("done")
    }

    ## All of the above can go into an SE.
    assays <- list(
        counts = counts, scaledCounts = counts / final_od, TPM = tpm,
        effectiveLength = eff_len
    )
    ## Handle a single sample case
    if (length(paths) == 1) {
        assays <- lapply(assays, as.matrix)
        counts <- as.matrix(counts)
    }
    rowDF <- DataFrame(
        length = lens[["Length"]], overdispersion = final_od, row.names = ids
    )
    colDF <- DataFrame(totals = colSums(counts), n_trans = n_trans)
    se <- SummarizedExperiment(assays = assays, rowData = rowDF, colData = colDF)
    metadata(se) <- list(resampleType = types)
    colnames(se) <- paths

    if (is(name_fun, "function")) colnames(se) <- name_fun(colnames(se))
    se

}

#' @importFrom matrixStats rowMeans2
#' @importFrom stats setNames median qf
#' @keywords internal
.overdispFromBoots <- function(paths, n_boot, n_trans, quants) {

    suf <- file.path("aux_info", "bootstrap", "bootstraps.gz")
    boot_files <- vapply(paths, file.path, character(1), suf)
    stopifnot(all(file.exists(boot_files)))

    boot_stats <- lapply(
        which(n_boot > 0), # There will always be at least one
        \(i){
            n_boot <- n_boot[[i]]
            n_trans <- n_trans[[i]]
            con <- gzcon(file(boot_files[[i]], open = "rb"))
            boots <- readBin(con, what = 'double', n = n_trans * n_boot)
            close(con)
            dim(boots) <- c(n_trans, n_boot)

            nm <- quants[[i]]$Name
            mn <- rowMeans2(boots)
            names(mn) <- nm

            df <- od <- setNames(rep_len(0, n_trans), nm)
            j <- mn > 0
            od[j] <- rowSums((boots[j,] - mn[j])^2) / mn[j]
            df[j] <- n_boot - 1L
            data.frame(transcript_id = nm, overdispersion = od, df = df)
        }
    )
    transcript_id <- c() ## R CMD check error avoidance
    boot_stats <- do.call("rbind", boot_stats)
    zeros <- subset(boot_stats, boot_stats$df == 0)
    boot_stats <- subset(boot_stats, boot_stats$df > 0)
    ## Unfortunately, this is much faster than splitting & lapplying
    overdispersion <- df <- NULL # R CMD check avoidance
    boot_sums <- dplyr::summarise(
        boot_stats,
        overdispersion = sum(overdispersion), df = sum(df), .by = transcript_id
    )
    boot_sums$od <- boot_sums$overdispersion / boot_sums$df
    df_med <- median(boot_sums$df)
    df_prior <- 3
    od_prior <- max(1, median(boot_sums$od) / qf(0.5, df_med, df_prior))
    od_num <- df_prior * od_prior + boot_sums$df * boot_sums$od
    od_denom <- df_prior + boot_sums$df
    boot_sums$od <- pmax(od_num / od_denom, 1)
    zero_ids <- setdiff(zeros$transcript_id, boot_sums$transcript_id)
    c(
        setNames(boot_sums$od, boot_sums$transcript_id),
        setNames(rep(od_prior, length(zero_ids)), zero_ids)
    )

}

.assayFromQuants <- function(x, var, fill = NA_real_) {

    df <- dplyr::bind_rows(x, .id = "sample")
    df <- as.data.frame(df[c("Name", "sample", var)])

    # wide <- tidyr::pivot_wider(
    #     df, names_from = "sample", values_from = var, values_fill = fill
    # )
    ## Try using reshape to avoid a tidyr dependency
    wide <- stats::reshape(
        df, idvar = "Name", direction = "wide", timevar = "sample",
        v.names = var
    )
    colnames(wide) <- gsub(paste0(var, "\\.", collapse = ""), "", colnames(wide))
    assay <- as.data.frame(lapply(wide[-1], \(x) {x[is.na(x)] <- fill; x}))
    rownames(assay) <- wide[[1]]
    as.matrix(assay)
}
