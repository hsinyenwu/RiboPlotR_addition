# One example for GTF files modification for RiboTaper

```
rm(list=ls())
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)

setwd("~/Desktop/atRTD3/")

#Import a gtf file
txdb <- makeTxDbFromGFF("~/Desktop/atRTD3/atRTD3_samples.gtf",format="gtf", dataSource="TAIR10",organism="Arabidopsis.thaliana")
gtf = import.gff("~/Desktop/atRTD3/atRTD3_samples.gtf")
grl = split(gtf, gtf$gene_id)  # grl a GRangesList
gr = unlist(grl) # Now gr is a GRanges object with multiple metacolumns such as seqnames, ranges, strand, type, phase, gene_id, transcript_id ... 

gr2 <- gr[gr$type=="exon"]
gr3 <- gr[gr$type=="CDS"]
gr2$gene_biotype <- ifelse(gr2$gene_id %in% gr3$gene_id, "protein_coding","ncRNA")
gr2$source="PBRI"
gr2$phase="."
gr2$score="."

gr3$gene_biotype <- "protein_coding"
gr3$score="."
gr3s <- split(gr3, gr3$transcript_id)

#Add phase to CDS rows
#https://support.bioconductor.org/p/101245/
addCdsPhase <- function(cds_by_tx) {
  cds_phase <- pc(rep(IntegerList(0), length(cds_by_tx)),
                  heads(cumsum(width(cds_by_tx)) %% 3L, n=-1L))
  unlisted_cds_by_tx <- unlist(cds_by_tx, use.names=FALSE)
  mcols(unlisted_cds_by_tx)$phase <- unlist(cds_phase, use.names=FALSE)
  relist(unlisted_cds_by_tx, cds_by_tx)
}

gr3s2 <- addCdsPhase(gr3s)
gr3s3 <- unlist(gr3s2)

gene_rows=range(grl)

gene_rows
s_gene_rows <- unlist(gene_rows)
s_gene_rows$type="gene"
s_gene_rows$source="PBRI"
s_gene_rows$phase="."
s_gene_rows$score="."
s_gene_rows$exon_number=NULL
s_gene_rows$gene_biotype <- ifelse(names(s_gene_rows) %in% gr3$gene_id, "protein_coding","ncRNA")
s_gene_rows$gene_id = names(s_gene_rows)

gr2$uniq_trans_id=NULL
gr2$exon_number=NULL
gr3s3$uniq_trans_id=NULL
gr3s3$exon_number=NULL

gtrC <- c(s_gene_rows,gr2,gr3s3)
gtrC2 <- split(gtrC,gtrC$gene_id)
gtrC2
length(gtrC2) #43

#remove the overlapping gene ranges
gtrC3 <- gtrC2[!grepl("-",names(gtrC2))]
length(gtrC3) #41

export(unlist(gtrC3), "~/Desktop/atRTD3/atRTD3_samples_modified.gtf",format="gff2")

At_gtf <- read.delim("~/Desktop/atRTD3/atRTD3_samples_modified.gtf",header=F,stringsAsFactors = F,quote = "",skip=3)
head(At_gtf)
At_gtf$V9 <- gsub("ID.*","",At_gtf$V9)
head(At_gtf,20)
write.table(At_gtf,file="~/Desktop/atRTD3/atRTD3_samples_modified.gtf", quote = F, col.names = F, row.names = F, sep = "\t",append = F)

# Change Chr1 to 1, ChrC to Pt, ChrM to Mt
# sed 's/Chr/chr/g' atRTD3_TS_21Feb22_transfix_wu.gtf > atRTD3_samples_modified.gtf
# sed 's/ChrC/Pt/g' atRTD3_TS_21Feb22_transfix_wu2.gtf > atRTD3_samples_modified.gtf
# sed 's/ChrM/Mt/g' atRTD3_TS_21Feb22_transfix_wu.gtf > atRTD3_samples_modified.gtf
# mv atRTD3_TS_21Feb22_transfix_wu2.gtf atRTD3_samples_modified.gtf

txdb <- makeTxDbFromGFF("~/Desktop/atRTD3/atRTD3_samples_modified.gtf",format="gtf", dataSource="TAIR10",organism="Arabidopsis thaliana")

```
