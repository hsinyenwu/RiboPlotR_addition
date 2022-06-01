# GTF files modification

```
rm(list=ls())
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
setwd("~/Desktop/atRTD3/")

#Import a gtf file
txdb <- makeTxDbFromGFF("~/Desktop/atRTD3/atRTD3_TS_21Feb22_transfix.gtf",format="gtf", dataSource="TAIR10",organism="Arabidopsis.thaliana")
gtf = import.gff("~/Desktop/atRTD3/atRTD3_TS_21Feb22_transfix.gtf")
grl = split(gtf, gtf$gene_id)  # grl a GRangesList
gr = unlist(grl) # Now gr is a GRanges object with multiple metacolumns such as seqnames, ranges, strand, type, phase, gene_id, transcript_id ... 

gr2 <- gr[gr$type=="gene"]
gr3 <- gr[gr$type=="exon"]
gr4 <- gr[gr$type=="CDS"]



```
