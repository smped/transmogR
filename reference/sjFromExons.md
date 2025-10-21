# Obtain Splice-Junctions from Exons and Transcripts

Using GRanges defining exons and transcripts, find the splice-junctions

## Usage

``` r
sjFromExons(
  x,
  rank_col = c("exon_number", "exon_rank"),
  tx_col = c("transcript_id", "tx_id"),
  extra_cols = "all",
  don_len = 8,
  acc_len = 5,
  as = c("GRanges", "GInteractions"),
  ...
)
```

## Arguments

- x:

  GRanges object with exons and transcripts. A column indicating the
  position (or rank) of each exon within the transcript must be
  included.

- rank_col:

  The column containing the position of each exons within the transcript

- tx_col:

  The column containing unique transcript-level identifiers

- extra_cols:

  Can be a vector of column names to return beyond rank_col and tx_col.
  By default all columns are returned (extra_cols = "all").

- don_len, acc_len:

  Length of donor and acceptor sites respectively

- as:

  Return as a set of GenomicRanges, or with each splice junction
  annotated as a GenomicInteraction

- ...:

  Not used

## Value

A GRanges object with requested columns, and an additional column,
'site', annotating each region as a donor or acceptor site.

Alternatively, by specifying as = "GInteractions", the junctions can be
returned with each splice junction annotated as a GenomicInteraction.
This can make the set of junctions easier to interpret for a given
transcript.

## Details

A canonical splice junction consists of a donor site and an acceptor
site at each end of an intron, with a branching site somewhere wthin the
intron. Canonical donor sites are 8nt long with the the first two bases
being exonic and the next 6 being derived form intronic sequences.
Canonical acceptor sites are 5nt long with the first four bases being
intronic and the final base being the first base of the next exon.

This functions uses each set of exons within a transcript to identify
both donor and acceptor sites. Branch sites are not identified.

## Examples

``` r
library(rtracklayer)
gtf_cols <- c(
  "transcript_id", "transcript_name", "gene_id", "gene_name", "exon_number"
)
gtf <- import.gff(
   system.file("extdata/gencode.v44.subset.gtf.gz", package = "transmogR"),
   feature.type = "exon", colnames = gtf_cols
)
sj <- sjFromExons(gtf)
sj
#> GRanges object with 730 ranges and 6 metadata columns:
#>         seqnames        ranges strand |        site     transcript_id
#>            <Rle>     <IRanges>  <Rle> | <character>       <character>
#>     [1]     chr1 182745-182752      + |       donor ENST00000624431.2
#>     [2]     chr1 183128-183132      + |    acceptor ENST00000624431.2
#>     [3]     chr1 183215-183222      + |       donor ENST00000624431.2
#>     [4]     chr1 183490-183494      + |    acceptor ENST00000624431.2
#>     [5]     chr1 183570-183577      + |       donor ENST00000624431.2
#>     ...      ...           ...    ... .         ...               ...
#>   [726]     chr1 810061-810068      - |       donor ENST00000590817.3
#>   [727]     chr1 810170-810174      - |    acceptor ENST00000447500.4
#>   [728]     chr1 817367-817374      - |       donor ENST00000447500.4
#>   [729]     chr1 827664-827671      - |       donor ENST00000635509.2
#>   [730]     chr1 827664-827671      - |       donor ENST00000634337.2
#>         transcript_name            gene_id       gene_name exon_number
#>             <character>        <character>     <character>   <integer>
#>     [1]    DDX11L17-201  ENSG00000279928.2        DDX11L17           1
#>     [2]    DDX11L17-201  ENSG00000279928.2        DDX11L17           2
#>     [3]    DDX11L17-201  ENSG00000279928.2        DDX11L17           2
#>     [4]    DDX11L17-201  ENSG00000279928.2        DDX11L17           3
#>     [5]    DDX11L17-201  ENSG00000279928.2        DDX11L17           3
#>     ...             ...                ...             ...         ...
#>   [726] ENST00000590817  ENSG00000230092.8 ENSG00000230092           1
#>   [727] ENST00000447500  ENSG00000290784.1 ENSG00000290784           2
#>   [728] ENST00000447500  ENSG00000290784.1 ENSG00000290784           1
#>   [729] ENST00000635509 ENSG00000230021.10 ENSG00000230021           1
#>   [730] ENST00000634337 ENSG00000230021.10 ENSG00000230021           1
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

## Or to simplify shared splice junctions across multiple transcripts
library(extraChIPs, quietly = TRUE)
chopMC(sj)
#> GRanges object with 196 ranges and 6 metadata columns:
#>         seqnames        ranges strand |        site
#>            <Rle>     <IRanges>  <Rle> | <character>
#>     [1]     chr1 182745-182752      + |       donor
#>     [2]     chr1 183128-183132      + |    acceptor
#>     [3]     chr1 183215-183222      + |       donor
#>     [4]     chr1 183490-183494      + |    acceptor
#>     [5]     chr1 183570-183577      + |       donor
#>     ...      ...           ...    ... .         ...
#>   [192]     chr1 808623-808627      - |    acceptor
#>   [193]     chr1 810061-810068      - |       donor
#>   [194]     chr1 810170-810174      - |    acceptor
#>   [195]     chr1 817367-817374      - |       donor
#>   [196]     chr1 827664-827671      - |       donor
#>                               transcript_id                 transcript_name
#>                             <CharacterList>                 <CharacterList>
#>     [1]                   ENST00000624431.2                    DDX11L17-201
#>     [2]                   ENST00000624431.2                    DDX11L17-201
#>     [3]                   ENST00000624431.2                    DDX11L17-201
#>     [4]                   ENST00000624431.2                    DDX11L17-201
#>     [5]                   ENST00000624431.2                    DDX11L17-201
#>     ...                                 ...                             ...
#>   [192] ENST00000447500.4,ENST00000590817.3 ENST00000447500,ENST00000590817
#>   [193] ENST00000447500.4,ENST00000590817.3 ENST00000447500,ENST00000590817
#>   [194]                   ENST00000447500.4                 ENST00000447500
#>   [195]                   ENST00000447500.4                 ENST00000447500
#>   [196] ENST00000635509.2,ENST00000634337.2 ENST00000635509,ENST00000634337
#>                                       gene_id                       gene_name
#>                               <CharacterList>                 <CharacterList>
#>     [1]                     ENSG00000279928.2                        DDX11L17
#>     [2]                     ENSG00000279928.2                        DDX11L17
#>     [3]                     ENSG00000279928.2                        DDX11L17
#>     [4]                     ENSG00000279928.2                        DDX11L17
#>     [5]                     ENSG00000279928.2                        DDX11L17
#>     ...                                   ...                             ...
#>   [192]   ENSG00000290784.1,ENSG00000230092.8 ENSG00000290784,ENSG00000230092
#>   [193]   ENSG00000290784.1,ENSG00000230092.8 ENSG00000290784,ENSG00000230092
#>   [194]                     ENSG00000290784.1                 ENSG00000290784
#>   [195]                     ENSG00000290784.1                 ENSG00000290784
#>   [196] ENSG00000230021.10,ENSG00000230021.10 ENSG00000230021,ENSG00000230021
#>           exon_number
#>         <IntegerList>
#>     [1]             1
#>     [2]             2
#>     [3]             2
#>     [4]             3
#>     [5]             3
#>     ...           ...
#>   [192]           3,2
#>   [193]           2,1
#>   [194]             2
#>   [195]             1
#>   [196]           1,1
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

## Splice Junctions can also be returned as a GInteractions object with
## anchorOne as the donor & anchorTwo as the acceptor sites
sjFromExons(gtf, as = "GInteractions")
#> GInteractions object with 365 interactions and 5 metadata columns:
#>         seqnames1       ranges1     seqnames2       ranges2 |     transcript_id
#>             <Rle>     <IRanges>         <Rle>     <IRanges> |       <character>
#>     [1]      chr1 182745-182752 ---      chr1 183128-183132 | ENST00000624431.2
#>     [2]      chr1 183215-183222 ---      chr1 183490-183494 | ENST00000624431.2
#>     [3]      chr1 183570-183577 ---      chr1 183736-183740 | ENST00000624431.2
#>     [4]      chr1 183900-183907 ---      chr1 183977-183981 | ENST00000624431.2
#>     [5]      chr1 724359-724366 ---      chr1 724714-724718 | ENST00000440782.3
#>     ...       ...           ... ...       ...           ... .               ...
#>   [361]      chr1 817367-817374 ---      chr1 774280-774284 | ENST00000447500.4
#>   [362]      chr1 827664-827671 ---      chr1 801163-801167 | ENST00000635509.2
#>   [363]      chr1 827664-827671 ---      chr1 805891-805895 | ENST00000634337.2
#>   [364]      chr1 720026-720033 ---      chr1 805891-805895 | ENST00000414688.6
#>   [365]      chr1 711705-711712 ---      chr1 808623-808627 | ENST00000635509.2
#>         transcript_name            gene_id       gene_name        sj
#>             <character>        <character>     <character> <integer>
#>     [1]    DDX11L17-201  ENSG00000279928.2        DDX11L17         1
#>     [2]    DDX11L17-201  ENSG00000279928.2        DDX11L17         2
#>     [3]    DDX11L17-201  ENSG00000279928.2        DDX11L17         3
#>     [4]    DDX11L17-201  ENSG00000279928.2        DDX11L17         4
#>     [5]       CICP3-201  ENSG00000229376.3           CICP3         1
#>     ...             ...                ...             ...       ...
#>   [361] ENST00000447500  ENSG00000290784.1 ENSG00000290784         1
#>   [362] ENST00000635509 ENSG00000230021.10 ENSG00000230021         1
#>   [363] ENST00000634337 ENSG00000230021.10 ENSG00000230021         1
#>   [364] ENST00000414688 ENSG00000230021.10 ENSG00000230021         2
#>   [365] ENST00000635509 ENSG00000230021.10 ENSG00000230021         2
#>   -------
#>   regions: 196 ranges and 0 metadata columns
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

```
