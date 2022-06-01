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

gr2$gene_biotype <- ifelse(gr2$gene_id %in% gr4$gene_id, "protein_coding","ncRNA")

gr3$gene_biotype <- ifelse(gr3$gene_id %in% gr4$gene_id, "protein_coding","ncRNA")
gr3$source="PBRI"
gr3$phase="."
gr3$score="."
#gr3$exon_number=NULL
#gr3$uniq_trans_id=NULL

gr4$gene_biotype <- "protein_coding"
gr4$score="."
#gr4$exon_number=NULL
#gr4$uniq_trans_id=NULL
gr4s <- split(gr4, gr4$transcript_id)

#Add phase to CDS rows
#https://support.bioconductor.org/p/101245/
addCdsPhase <- function(cds_by_tx) {
  cds_phase <- pc(rep(IntegerList(0), length(cds_by_tx)),
                  heads(cumsum(width(cds_by_tx)) %% 3L, n=-1L))
  unlisted_cds_by_tx <- unlist(cds_by_tx, use.names=FALSE)
  mcols(unlisted_cds_by_tx)$phase <- unlist(cds_phase, use.names=FALSE)
  relist(unlisted_cds_by_tx, cds_by_tx)
}

gr4s2 <- addCdsPhase(gr4s)
gr4s3 <- unlist(gr4s2)



```
