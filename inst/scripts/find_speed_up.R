#' This is an attempt to revamp the transmogrify function for increased speed.
#' Instead of inserting SNPs then InDels in the reference, as relevant for
#' each transcript, try:
#'
#' 1. inserting into the reference, creating a map of old-new co-ordinates
#' 2. extract transcriptSeqs as per-usual
#'
#' - This will also produce a variant-modified genome
#' - Chunking by chromosome by also improve performance, but knowing which steps
#'   to chunk my be non-trivial
#'

library(tidyverse)
library(VariantAnnotation)
library(GenomicRanges)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(plyranges)
library(transmogR)
library(rtracklayer)
library(GenomicFeatures)
library(parallel)

## Load the relevant reference, then modify
ref <- BSgenome.Hsapiens.UCSC.hg38
chr <- paste0("chr", c(1:22, "X", "Y", "M"))
# chr <- paste0("chr", c(11:15))

## Start by loading a set of variants, subset to chr10 for testing
var <- read_rds("~/TKI/DBPAInT/data/rds/1000GP_SNV_INDEL_panhuman.rds") |>
    subset(seqnames %in% chr) |>
    transmogR:::.checkOverlapVars(ol_vars = "none")
new_ref <- genomogrify(ref, var, names = chr)

## Load the GTF
gtf <- read_rds("~/TKI/DBPAInT/data/rds/gencode.v44.rds") %>%
    subset(grepl("chr", seqnames)) %>%
    splitAsList(.$type)
exons_by_trans <- gtf$exon %>% splitAsList(.$transcript_id)
trans_seq <- extractTranscriptSeqs(ref, exons_by_trans)

## Test the new method
new_exon <- shiftByVar(gtf$exon, var, mc.cores = 2)
new_exons_by_trans <- splitAsList(new_exon, new_exon$transcript_id)
new_trans_seq <- extractTranscriptSeqs(new_ref, new_exons_by_trans)
mean(new_trans_seq == trans_seq)
which(new_trans_seq != trans_seq) |> head()


# Now Figure it out by changing co-ords
change <- Rle(0L, nchar(ref))
change[start(var_map)] <- var_map$change
shift <- cumsum(change)
new_exon <- GRanges(
    seqnames = seqnames(gtf$exon),
    IRanges(
        start = start(gtf$exon) + shift[start(gtf$exon)],
        end = end(gtf$exon) + shift[end(gtf$exon)] ,
    ),
    strand = strand(gtf$exon),
    seqinfo = seqinfo(gtf$exon)
)
mcols(new_exon) <- mcols(gtf$exon)
new_exons_by_trans <- new_exon %>% splitAsList(.$transcript_id)
new_trans_seq <- extractTranscriptSeqs(new_ref, new_exons_by_trans)
trans_seq
new_trans_seq
sum(new_trans_seq == trans_seq)
sum(new_trans_seq != trans_seq)

## Now everything works. Benchmark...
## The new approach
peakRAM::peakRAM(
    # {
        ## This is now the complete process & it's 8-10x faster
        ## Adding tags has not yet been incorporated though and that will add
        ## some time to the process
        new_ref <- genomogrify(ref, var),
        ## Create the map of new-old
        split_var <- split(var, varTypes(var)),
        split_var <- as.list(split_var),
        split_var$SNV <- NULL,
        split_var$Deletion <- GPos(split_var$Deletion),
        split_var$Deletion$id <- subjectHits(findOverlaps(split_var$Deletion, var)),
        split_var$Deletion$change <- -1 * c(0, diff(start(split_var$Deletion))),
        ## This handles directly neighbouring deletions & also removes the first position
        ## of any deletion
        keep <- c(FALSE, diff(split_var$Deletion$id) == 0),
        split_var$Deletion <- split_var$Deletion[keep],
        split_var$Insertion$change <- nchar(split_var$Insertion$ALT) - nchar(split_var$Insertion$REF),
        map <- GRangesList(lapply(split_var, GRanges)),
        map <- sort(unlist(map)),
        change <- rep_len(0, nchar(ref)),
        change[start(map)] <- map$change,
        shift <- cumsum(change),
        new_exon <- GRanges(
            seqnames = seqnames(gtf$exon),
            IRanges(
                start = start(gtf$exon) + shift[start(gtf$exon)],
                end = end(gtf$exon) + shift[end(gtf$exon)] ,
            ),
            strand = strand(gtf$exon),
            seqinfo = seqinfo(gtf$exon)
        ),
        mcols(new_exon) <- mcols(gtf$exon),
        new_exons_by_trans <- new_exon %>% splitAsList(.$transcript_id),
        new_trans_seq <- extractTranscriptSeqs(new_ref, new_exons_by_trans)
    # }
) -> pr_df2
# user  system elapsed
# 8.607   0.363   8.990
## The old approach
system.time(
    {
        old_trans_seq <- transmogrify(ref, var, unlist(exons_by_trans))
        new_ref <- genomogrify(ref, var)
    }
)
#   user  system elapsed
# 75.122   0.095  75.336
sum(old_trans_seq != new_trans_seq)
# [1] 0
## OK, now to put it into a function
## Should the function take a modified reference, or modify it
## If modifying, how should it be returned in the final object?
## RAM may be an issue
## Maybe modifying internally & then running `genomogrify()` separately is
## still best...?
