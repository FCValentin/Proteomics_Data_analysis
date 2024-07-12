### READ FILES ###

# Import some home functions
source("https://gitlab.univ-nantes.fr/E114424Z/veneR/raw/master/loadFun.R?inline=false")

# Set working directory

#!/usr/bin/env Rscript
folder <- commandArgs(trailingOnly=TRUE)
setwd(folder)

# Import files to use for the analyze
filenames <- list.files(path="Comparaison", pattern="*.tsv", full.names=TRUE, all.files = TRUE, recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE)
print("Files used for the analysis :")
print(filenames)

# Import proteins background
ProtNames <- row.names(lire("data/dataMatrix.tsv"))
print(paste("Total of proteins :", length(ProtNames)))

# Create data matrix of protein background
InitTable<-data.frame(row.names = ProtNames)
RawComparisonTable <- InitTable[order(row.names(InitTable)),]

print("Matrix initialization...")

# Initialization of the table with 0 everywhere
for (i in 1:length(filenames)){
  RawComparisonTable<-cbind(RawComparisonTable,0)
  colnames(RawComparisonTable)[i]<-strsplit(strsplit(filenames[i],split = "/")[[1]][2],split = ".t")[[1]][1]
}

print("Matrix formating...")

# Replace 0 values by AVG_Log2Ratio value for DE Genes in each comparison
for (i in 1:length(filenames)){
  file<-read.table(filenames[i])
  colnames(file)<-c("Protein","AVG_Log2_Ratio","Absolute_AVG_Log2_Ratio","PValue","QValue","%ofRatios","UniProtIds","Genes","ProteinName","%UniqueTotalPeptide","%Change","Ratio")
  ordered_file <- file[order(row.names(file)),]
  RawComparisonTable[row.names(ordered_file),i]<-ordered_file$AVG_Log2_Ratio
}

print("Output the matrix in results Folder as ComparisonTable.tsv ...")

ecrire(RawComparisonTable,"results/ComparisonTable.tsv")

BinaryComparisonTable <- RawComparisonTable
BinaryComparisonTable[BinaryComparisonTable>0] <- 1
BinaryComparisonTable[BinaryComparisonTable<0] <- -1

print("Output the binary matrix in results Folder as BinaryComparisonTable.tsv ...")

ecrire(BinaryComparisonTable,"results/BinaryComparisonTable.tsv")

