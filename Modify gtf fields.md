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

gtr <- c(gr2,gr3,gr4s3)
gtr2 = split(gtr, gtr$gene_id)
length(gtr2) #40929
gtr2l <- unlist(gtr2)
export(gtr2l, "~/Desktop/atRTD3/atRTD3_TS_21Feb22_transfix_wu_test.gtf",format="gff2")

At_gtf <- read.delim("~/Desktop/atRTD3/atRTD3_TS_21Feb22_transfix_wu_test.gtf",header=F,stringsAsFactors = F,quote = "",skip=3)
head(At_gtf)

At_gtf$V9 <- gsub("ID.*","",At_gtf$V9)

write.table(At_gtf,file="~/Desktop/atRTD3/atRTD3_TS_21Feb22_transfix_wu_test.gtf", quote = F, col.names = F, row.names = F, sep = "\t",append = F)
txdb2 <- makeTxDbFromGFF("~/Desktop/atRTD3/atRTD3_TS_21Feb22_transfix_wu_test.gtf",format="gtf", dataSource="TAIR10",organism="Arabidopsis.thaliana")

# gene_rows <- lapply(seq_along(grl), function(x) unlist(range(grl[x])))

A=range(grl)
#length(grl) [1] 40929
for(x in 1:40929) A[[x]]$gene_id=names(grl)[x]
gene_rows=A
save(gene_rows,file ="~/Desktop/atRTD3/gene_rows.RData")
load("~/Desktop/atRTD3/gene_rows.RData")

gene_rows
s_gene_rows <- unlist(gene_rows)
s_gene_rows$type="gene"
s_gene_rows$source="PBRI"
s_gene_rows$phase="."
s_gene_rows$score="."
s_gene_rows$exon_number=NULL
s_gene_rows$gene_biotype <- ifelse(s_gene_rows$gene_id %in% gr4$gene_id, "protein_coding","ncRNA")

gtrC <- c(s_gene_rows,gr3,gr4s3)
gtrC2 <- split(gtrC,gtrC$gene_id)
gtrC2
length(gtrC2) #40929
#remove the overlapping gene ranges
gtrC3 <- gtrC2[!grepl("-",names(gtrC2))]
length(gtrC3) #39732

export(unlist(gtrC3), "~/Desktop/atRTD3/atRTD3_TS_21Feb22_transfix_wu.gtf",format="gff2")

At_gtf <- read.delim("~/Desktop/atRTD3/atRTD3_TS_21Feb22_transfix_wu.gtf",header=F,stringsAsFactors = F,quote = "",skip=3)
head(At_gtf)
At_gtf$V9 <- gsub("ID.*","",At_gtf$V9)
head(At_gtf,20)
write.table(At_gtf,file="~/Desktop/atRTD3/atRTD3_TS_21Feb22_transfix_wu.gtf", quote = F, col.names = F, row.names = F, sep = "\t",append = F)

# Change Chr1 to 1, ChrC to Pt, ChrM to Mt
# sed 's/Chr/chr/g' atRTD3_TS_21Feb22_transfix_wu.gtf > atRTD3_TS_21Feb22_transfix_wu2.gtf
# sed 's/ChrC/Pt/g' atRTD3_TS_21Feb22_transfix_wu2.gtf > atRTD3_TS_21Feb22_transfix_wu.gtf
# sed 's/ChrM/Mt/g' atRTD3_TS_21Feb22_transfix_wu.gtf > atRTD3_TS_21Feb22_transfix_wu2.gtf
# mv atRTD3_TS_21Feb22_transfix_wu2.gtf atRTD3_TS_21Feb22_transfix_wu.gtf

txdb <- makeTxDbFromGFF("~/Desktop/atRTD3/atRTD3_TS_21Feb22_transfix_wu.gtf",format="gtf", dataSource="TAIR10",organism="Arabidopsis thaliana")

```
