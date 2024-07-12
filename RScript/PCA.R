### PARAMETERS TO CHANGE ###
directory<-"~/These/Analyses/NolwennProject" #Enter the path to the directory to use
Matrix<-"data/dataMatrix_wash_groups_newSettings.tsv" # Enter the folder and path to the matrix data file that you want to read, from your working directory position
samples<-"data/SampleAnnot.tsv" # Enter the folder and path to the samples data file that you want to read, from your working directory position
NbPCA_Axis<-4 # Enter the max number of PCA axis you want to show

# Import some (home) functions
source("https://gitlab.univ-nantes.fr/E114424Z/veneR/raw/master/loadFun.R?inline=false")
library(ggplot2)

# Set working directory
setwd(directory)

# Import proteins background
data <- lire(Matrix)
sampleAnnot <- lire(samples)

# Check if there are problems on sample names
sampleExprInAnnot<-rn(sampleAnnot)[!rn(sampleAnnot)%in%cn(ordered_data)]; if(len(sampleExprInAnnot)>0){warning("Error, these samples don't exist in expression data, please check sample names"); print(sampleExprInAnnot)}
sampleAnnotInExpr<-cn(ordered_data)[!cn(ordered_data)%in%rn(sampleAnnot)]; if(len(sampleAnnotInExpr)>0){warning("these samples don't exist in sample table, they will be removed from expression data"); print(sampleAnnotInExpr)}
ordered_data <- data[order(row.names(data)),row.names(sampleAnnot)]
print("Sample names match is ok")

##### PCA Quality control#####
print(paste("Total of clusters :",NbPCA_Axis))

Batch <- "Batch"
Feader <- "Feader"
Replication <- "Geminin"

acp<-ACP(ordered_data[,which(sampleAnnot$Tissue=="Sperm")]) 
pdf(file = paste0("Figures/PCA_Sperm",NbPCA_Axis,"_axis.pdf"),width=10,height=10)
barplot(acp$percentVar,names.arg = round(acp$percentVar*100,2),main = "Contribution of each componant in PCA")

print("Performing PCA....")

for(i in 1:(NbPCA_Axis-1)){
  for(j in (i+1):NbPCA_Axis){
    acp2d(acp,group = as.factor(sampleAnnot[which(sampleAnnot$Tissue=="Sperm"),Batch]),plotText = TRUE,pointSize = 4,comp = c(i,j),main="PCA Batch",fixedCoord = F)
    acp2d(acp,group = as.factor(sampleAnnot[which(sampleAnnot$Tissue=="Sperm"),Batch]),pointSize = 2,comp = c(i,j),main="PCA Batch",fixedCoord = F)
    
    acp2d(acp,group = sampleAnnot[which(sampleAnnot$Tissue=="Sperm"),Feader],plotText = TRUE,pointSize = 4,comp = c(i,j),main="PCA Feader",fixedCoord = F)
    acp2d(acp,group = sampleAnnot[which(sampleAnnot$Tissue=="Sperm"),Feader],pointSize = 2,comp = c(i,j),main="PCA Feader",fixedCoord = F)
    
    acp2d(acp,group = sampleAnnot[which(sampleAnnot$Tissue=="Sperm"),Replication],plotText = TRUE,pointSize = 4,comp = c(i,j),main="PCA Replication",fixedCoord = F)
    acp2d(acp,group = sampleAnnot[which(sampleAnnot$Tissue=="Sperm"),Replication],pointSize = 2,comp = c(i,j),main="PCA Replication",fixedCoord = F)
    
    acp2d(acp,pointSize = 2,comp = c(i,j),plotVars = TRUE, plotText = TRUE,fixedCoord = F)
  }
}
dev.off()
print("End of the script, output in Figures folder")



