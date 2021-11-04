## Here we describe how to generate gtf files for uORFs

Remember there could be multiple uORFs for one 5'UTR, and there will be multiple gtf files generate in that situation.

### Generate gtf files for uORFs using genomic sequences and annotation for 5'UTRs and transcripts   

Note:
1. Here we only consider uORFs start with AUG but you can modify the start sequences in the findUORFs (see ORFik package for detail) function
2. Requirements: (1) a fasta file for genome sequence (2) a gtf files with CDS information
3. The following code works in Mac OS and Linux, but might require slight modification for Windows system
4. The compressed files for: TAIR10_chr_1.fa and Araport11+CTRL_20181206.gtf are attached so you can try this function out.

```
# Make uORF GTFs using genomic sequences and annotation for 5'UTRs and transcripts for one transcript
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
  #Find all possible uORFs
  uORF <- findUORFs(fiveUTRs=fiveUTR[y], startCodon="ATG", fa=FA, longestORF=T, cds=CDS[y], restrictUpstreamToTx=T)
  #Export the exon info for the uORF containing transcript to a temp gtf file
  export.gff(exonRanges[y],con=paste0(z, "uORF_exon_ranges.gtf"), format="gff2")
  #Read in the exon info. Now it is in the gtf format
  EXONu <- read.delim(paste0(z,"uORF_exon_ranges.gtf"), skip=4, header=F)
  #Modify the column 3 and 9
  EXONu$V3 <- "mRNA"
  EXONu$V9 <- paste0("gene_id ","\"",x,"\"; ", "transcript_id ","\"",y,"\";")
  
  for (i in seq_along(uORF)) {
  #Export the uORF CDS range for uORF i (i is the number of a specific uORF)
    export.gff(uORF[i], con=paste0(z,"uORF_CDS_ranges.gtf"), format="gff2")
    #Read in the CDS range
    CDSu <- read.delim(paste0(z,"uORF_CDS_ranges.gtf"), skip=4, header=F)
    #Modify column 3 and 9
    CDSu$V3 <- "CDS"
    CDSu$V9 <- paste0("gene_id ", "\"",x, "\"; ", "transcript_id ","\"",y,"\";")
    rbind(EXONu,CDSu)
    write.table(rbind(EXONu,CDSu), paste0(z,y,".",i,".gtf"), quote = F, col.names = F, row.names = F, sep = "\t",append = F)
  }
}

#Example:
generate_uORFgtfs(x="AT1G01060",y="AT1G01060.1",z="~/Desktop/test/") #Should generate five uORF gtf files.

```



Citation: Please cite our paper: https://www.biorxiv.org/content/10.1101/694646v1 if you use information provided here for your research. Thanks!
