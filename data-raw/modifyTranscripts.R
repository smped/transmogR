#' First attempt at modifying transcripts using the 1000GP variants

library(tidyverse)
library(rtracklayer)
library(plyranges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BiocParallel)
source("code/putSNPsIntoSeq.R")
#' There is a function `injectSNPs` which can be used easily here
#' However, the SNPs must be a SNPlocs object and these currently appear to
#' only be available as packages
hg38 <- BSgenome.Hsapiens.UCSC.hg38
sq <- seqinfo(hg38)
gtf <- here::here("data", "annotations", "gencode.v44.primary_assembly.annotation.gtf.gz") %>%
  import.gff() %>%
  subset(seqnames %in% seqlevels(sq)) %>%
  keepSeqlevels(seqlevelsInUse(.)) %>%
  split(.$type)
seqinfo(gtf) <- sq


variant_pos <- here::here("data", "annotations", "snp_indel_coords.txt.gz") %>%
  read_tsv(
    col_names = c("CHROM", "POS", "ID", "REF", "ALT"), col_types = "ciccc",
    comment = "#"
  ) %>%
  dplyr::filter(!str_detect(ID, "HGSV")) %>%
  mutate(
    ref_width = str_count(REF), alt_width = str_count(ALT),
    type = case_when(
      ref_width > alt_width ~ "Deletion",
      ref_width < alt_width ~ "Insertion",
      ref_width + alt_width == 2 ~ "SNV"
    ),
    width = ref_width,
    end = POS + width - 1,
    inframe = abs(ref_width - alt_width) %% 3 == 0
  )
gr <- variant_pos %>%
  makeGRangesFromDataFrame(
    keep.extra.columns = TRUE,
    seqnames.field = "CHROM", start.field = "POS", end.field = "end",
    ignore.strand = TRUE#,
    # seqinfo = sq
  )
rm(variant_pos)

#' Try setting up the entire genome as a DNAStringSet, then extracting transcripts
hg38_mod <- putSNPsIntoSeq(hg38, subset(gr, type == "SNV"), "ALT")

#" Extract the SNP-modified transcripts
trByExon <- splitAsList(gtf$exon, gtf$exon$transcript_id)
trSeqByExon <- getSeq(hg38_mod, trByExon)
trSeq <- endoapply(trSeqByExon[1:5000], unlist)
trWithIns <- gtf$exon %>%
  subset(transcript_id %in% names(trSeq)) %>%
  select(transcript_id, exon_number) %>%
  split(.$transcript_id) %>%
  subsetByOverlaps(subset(gr, type != "SNV"))



