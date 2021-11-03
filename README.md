### Here we describe how to generate gtf files for uORFs

Remember there could be multiple uORFs for one 5'UTR so there will be multiple gtf files generate:

'''
# Make uORF GTFs for one transcript
# rm(list=ls())
library(GenomicRanges)
library(GenomicFeatures)
library(ORFik)
library(rtracklayer)

# Load the genome sequence file (e.g. TAIR10 genome sequence):
FA <- FaFile("~/Desktop/Leaky_scanning/TAIR10_chr_all_2.fas")
# Make the TxDb file (Use a gtf from Araport11)
txdb <- makeTxDbFromGFF("~/Desktop/CTRL_v1/Araport11+CTRL_20181206.gtf",format="gtf", dataSource="Araport11",organism="Arabidopsis")
# Extract 5'UTR ranges
fiveUTR <- fiveUTRsByTranscript(txdb, use.names=T)
# Extract CDS ranges
CDS <- cdsBy(txdb,by="tx", use.names=T)
# Extract exon ranges
exonRanges <- exonsBy(txdb,by="tx", use.names=T)

# generate_uORFgtf, parameter x is the transcript name
generate_uORFgtfs <-function(x,y,z){
  uORF <<- findUORFs(fiveUTRs=fiveUTR[y],fa=FA,longestORF=T,cds=CDS[y],restrictUpstreamToTx=T)
  export.gff(exonRanges[y],con=paste0(z,"uORF_exon_ranges.gtf"),format="gff2")
  EXONu <- read.delim(paste0(z,"uORF_exon_ranges.gtf"),skip=4,header=F)
  EXONu$V3 <- "mRNA"
  EXONu$V9 <- paste0("gene_id ","\"",x,"\"; ", "transcript_id ","\"",y,"\";")
  for (i in seq_along(uORF)) {
    export.gff(uORF[i],con=paste0(z,"uORF_CDS_ranges.gtf"),format="gff2")
    CDSu <- read.delim(paste0(z,"uORF_CDS_ranges.gtf"),skip=4,header=F)
    CDSu$V3 <- "CDS"
    CDSu$V9 <- paste0("gene_id ","\"",x,"\"; ", "transcript_id ","\"",y,"\";")
    rbind(EXONu,CDSu)
    write.table(rbind(EXONu,CDSu),paste0(z,y,".",i,".gtf"),quote = F, col.names = F, row.names = F, sep = "\t",append = F)
  }
}

generate_uORFgtfs("AT1G01060","AT1G01060.1","~/Desktop/test/")
'''

Citation: Please cite our paper: https://www.biorxiv.org/content/10.1101/694646v1 if you use information provided here for your research. 
