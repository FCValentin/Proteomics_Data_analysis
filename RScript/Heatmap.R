### PARAMETERS TO CHANGE ###
directory<-"~/These/Analyses/NolwennProject" #Enter the path to the directory to use
Matrix<-"data/dataMatrix_wash_groups_newSettings.tsv" # Enter the folder and path to the matrix data file that you want to read, from your working directory position
samples<-"data/SampleAnnot.tsv" # Enter the folder and path to the samples data file that you want to read, from your working directory position
NbPCA_Axis<-4 # Enter the max number of PCA axis you want to show
ListOfProt <- "ListOfDemandes/list_orthologs_EPI_and_TF_SEGM_DE/list_orthologs_EPI_and_TF_SEGM_DE.txt" # Enter the folder and path to the Protein list data file to plot, from your working directory position
filenames <- "Comparaison/Sperm_prog.tsv" # Enter the folder and path to the compairaison list data file to have genename, from your working directory position

# 1 = Prog / Sperm / proggem- / spermgem- 
# 2 = Prog / spermgem- / proggem- / spermgem+ / proggem+
# 3 = proggem- / spermgem- / proggem+ / spermgem+ 
# 4 = sperm / sperm DE gem-
# 5 = sperm / sperm DE gem- / sperm DE gem+
# 6 = prog / prog DE gem-
# 7 = prog / prog DE gem- / prog DE gem+
# 8 = Prog / Sperm
# 9 = personal ordering
WhatToCompare<-6

# Import some home functions
source("https://gitlab.univ-nantes.fr/E114424Z/veneR/raw/master/loadFun.R?inline=false")
library(circlize)
library(ComplexHeatmap)
library(Matrix)

# Set working directory
setwd(directory)

# Import proteins background
data <- lire(Matrix)
sampleAnnot <- lire(samples)

# Check if there are problems on sample names
sampleExprInAnnot<-rn(sampleAnnot)[!rn(sampleAnnot)%in%cn(data)]; if(len(sampleExprInAnnot)>0){warning("Error, these samples don't exist in expression data, please check sample names"); print(sampleExprInAnnot)}
sampleAnnotInExpr<-cn(data)[!cn(data)%in%rn(sampleAnnot)]; if(len(sampleAnnotInExpr)>0){warning("these samples don't exist in sample table, they will be removed from expression data"); print(sampleAnnotInExpr)}
ordered_data <- data[order(row.names(data)),row.names(sampleAnnot)]
print("Sample names match is ok")

#### Script for candidates Genes list ####
DEProt<-read.table(ListOfProt)
ProtBackground<-unique(DEProt$V1)

if(WhatToCompare==1){
  LinesToCompare<-c(1,2,3,4,21,22,23,24,5,6,7,8,25,26,27,28)   
}else if(WhatToCompare==2){
  LinesToCompare<-c(1,2,3,4,21,22,23,24,5,6,7,8,25,26,27,28,9,10,11,12,29,30,31,32)   
}else if(WhatToCompare==3){
  LinesToCompare<-c(5,6,7,8,25,26,27,28,9,10,11,12,29,30,31,32) 
}else if(WhatToCompare==4){
  LinesToCompare<-c(21,22,23,24,25,26,27,28)   
}else if(WhatToCompare==5){
  LinesToCompare<-c(21,22,23,24,25,26,27,28,29,30,31,32)  
}else if(WhatToCompare==6){
  LinesToCompare<-c(1,2,3,4,5,6,7,8)  
}else if(WhatToCompare==7){
  LinesToCompare<-c(1,2,3,4,5,6,7,8,9,10,11,12) 
}else if(WhatToCompare==8){
  LinesToCompare<-c(1,2,3,4,21,22,23,24)
}else{
  LinesToCompare<-WhatToCompare
}

Data<-ordered_data[which(row.names(ordered_data)%in%ProtBackground),LinesToCompare]
file<-lire(filenames)
mat<-rowScale(Data,center = TRUE,scaled = TRUE)
hclustGeneDE<-unsupervisedClustering(mat,transpose = F,nboot=30,bootstrap = FALSE,method.dist="pearson")
quantile.expr<-quantile(unlist(mat),seq(0,1,.01))
colHA<-colorRamp2(c(quantile.expr[2],0,quantile.expr[100]),c("blue","white","red"))
Ht<-Heatmap(mat,row_labels = file[row.names(mat),"Genes"], row_names_gp = autoGparFontSizeMatrix(nrow(mat)),
            cluster_rows = hclustGeneDE,col = colHA,
            cluster_columns = FALSE,name=paste0("(",nrow(mat)," Proteins)"),column_names_gp = autoGparFontSizeMatrix(ncol(mat)) )

max<-nrow(mat)
HeigthScaling<-10
if(max<400){
  HeigthScaling<-10;
}else if(max<1000){
  HeigthScaling<-15;
}else if(max<2000){
  HeigthScaling<-20;
}else if(max<4000){
  HeigthScaling<-40;
}else if(max<6000){
  HeigthScaling<-50;
}else{
  HeigthScaling<-60;
}

pdf(paste0("Figures/",paste0("(",nrow(mat)," Proteins)"),".pdf"),width = 10,height=HeigthScaling)
print(Ht)
dev.off()
ecrire(mat,paste0("results/",paste0("(",nrow(mat)," Proteins)"),".tsv"))
