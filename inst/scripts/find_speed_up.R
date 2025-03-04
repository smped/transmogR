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

## Start by loading a set of variants, subset to chr10 for testing
var <- read_rds("~/TKI/DBPAInT/data/rds/1000GP_SNV_INDEL_panhuman.rds") |>
    subset(seqnames == "chr10")

## Load the relevant reference, then subset to chr10 & modify
ref <- getSeq(BSgenome.Hsapiens.UCSC.hg38, names = "chr10") |>
    as("DNAStringSet") |>
    setNames("chr10")
new_ref <- genomogrify(ref, var)

## Create the map of new-old
var_map <- var %>%
    subset(nchar(REF) != nchar(ALT)) %>%
    mutate(
        change = nchar(ALT) - nchar(REF),
        cumsum_change = cumsum(change),
        new_start = start + c(0, cumsum_change[-length(.)]),
        new_end = end + cumsum_change
    )

## Load the GTF
gtf <- read_rds("~/TKI/DBPAInT/data/rds/gencode.v44.rds") %>%
    subset(seqnames == "chr10") %>%
    splitAsList(.$type)
exons_by_trans <- gtf$exon %>%
    splitAsList(.$transcript_id)

## extractTranscriptSeqs
trans_seq <- extractTranscriptSeqs(ref, exons_by_trans)

## Now Figure it out by changing co-ords
change <- rep_len(0L, nchar(ref))
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
sum(new_trans_seq == trans_seq) # [1] 57 ## Nope
## ENST00000016171.6 is off by 9
## ENST00000020673.6 is irreconcilable
## ENST00000173785.4 is off by -1

tail(shift)
# [1] -968 -968 -968 -968 -968 -968
nchar(ref) - nchar(new_ref)
# [1] 956
## So there's an error of ~142t using this strategy which is also present in the map!!
## Or is genomogrify wrong? I've checked the transcripts with transmogrify...

## Maybe handling each position instead of each range
var %>%
    subset(nchar(REF) != nchar(ALT)) %>%
    mutate(
        shift = nchar(ALT) - nchar(REF)
    )

## The issues are problems caused by the following
var[c(88226, 88227)]
# GRanges object with 2 ranges and 2 metadata columns:
#                        seqnames              ranges strand |         REF         ALT
#                           <Rle>           <IRanges>  <Rle> | <character> <character>
#     10:114468420:GCC:G    chr10 114468420-114468422      * |         GCC           G
#    10:114468422:C:CTAT    chr10           114468422      * |           C        CTAT
# -------
#     seqinfo: 24 sequences from an unspecified genome
## How do we delete the GCC (replacing with G), then include the C from the CTAT in the insertion
## Realistically, the GCC is being replaced by GTAT



var[c(2708, 2709)]
# GRanges object with 2 ranges and 2 metadata columns:
#                             seqnames          ranges strand |         REF         ALT
#                                <Rle>       <IRanges>  <Rle> | <character> <character>
#     10:2569741:TTGTGTACTC:T    chr10 2569741-2569750      * |  TTGTGTACTC           T
#           10:2569750:C:CAAA    chr10         2569750      * |           C        CAAA
# -------
#     seqinfo: 24 sequences from an unspecified genome
## Similarly, the TTGTGTACTC is being replaced by TAAA
