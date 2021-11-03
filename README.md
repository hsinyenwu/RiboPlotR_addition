### Here we describe how to generate gtf files for uORFs

Remember there could be multiple uORFs for one 5'UTR, and there will be multiple gtf files generate in that situation.

Note:
1. Here we only consider uORFs start with AUG but you can modify the start sequences in the findUORFs (see ORFik package) function
2. Requirements: (1) a fasta file for genome sequence (2) a gtf files with CDS information
3. The following code works in Mac OS and Linux, but might require slight modification for Windows system
4. The following files TAIR10_chr_1.fa and Araport11+CTRL_20181206.gtf are attached.

```
# Make uORF GTFs for one transcript
# rm(list=ls())
library(GenomicRanges)
library(GenomicFeatures)
library(Rsamtools)
library(ORFik)
library(rtracklayer)

# Load the genome sequence file (e.g. TAIR10 genome sequence):
FA <- FaFile("~/Desktop/TAIR10_chr_1.fa")
# Make the TxDb file (Use a gtf from Araport11)
txdb <- makeTxDbFromGFF("~/Desktop/Araport11+CTRL_20181206.gtf",format="gtf", dataSource="Araport11",organism="Arabidopsis")
# Extract 5'UTR ranges
fiveUTR <- fiveUTRsByTranscript(txdb, use.names=T)
# Extract CDS ranges
CDS <- cdsBy(txdb,by="tx", use.names=T)
# Extract exon ranges
exonRanges <- exonsBy(txdb,by="tx", use.names=T)

# generate_uORFgtfs function
# x is the gene name
# y is the transcript name
# z is the path for the new gtf files

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

generate_uORFgtfs("AT1G01060","AT1G01060.1","~/Desktop/test/") #Should generate five uORF gtf files
```

Citation: Please cite our paper: https://www.biorxiv.org/content/10.1101/694646v1 if you use information provided here for your research. Thanks!
