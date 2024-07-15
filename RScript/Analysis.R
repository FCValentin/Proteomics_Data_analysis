### READ FILES ###

# Import some home functions
source("https://gitlab.univ-nantes.fr/E114424Z/veneR/raw/master/loadFun.R?inline=false")

library(circlize)
library(ComplexHeatmap)
library(Matrix)

Batch <- "Batch"
Feader <- "Feader"
Replication <- "Geminin"
CellType <- "Tissue"

# Set working directory
#!/usr/bin/env Rscript
folder <- commandArgs(trailingOnly=TRUE)[1]
setwd(folder)
setwd("~/These/Analyses/NolwennProject")

# Import proteins background
data <- lire("data/dataMatrix_wash_groups_newSettings.tsv")
sampleAnnot <- lire("data/SampleAnnot.tsv")
ordered_data <- data[order(row.names(data)),row.names(sampleAnnot)]


if((length(commandArgs(trailingOnly=TRUE))-1)>1){
  Variables<-list()
  for (j in 1:(length(commandArgs(trailingOnly=TRUE))-1)){
    Variables[[j]] <- commandArgs(trailingOnly=TRUE)[j+1]
    if(!Variables[[j]]%in%colnames(sampleAnnot)){print(Variables[[j]]); stop("Error, previous parameter is not named in the sampleAnnot file, please check this condition")}
  }
}
print("Conditions to show on heatmap are :")
print(unlist(Variables))

# Import files to use for the analyze
filenames <- list.files(path="Comparaison", pattern="*.tsv", full.names=TRUE, all.files = TRUE, recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE)
print("Files used for the analysis :")
print(filenames)

# Check if there are problems on sample names
sampleExprInAnnot<-rn(sampleAnnot)[!rn(sampleAnnot)%in%cn(ordered_data)]; if(len(sampleExprInAnnot)>0){warning("Error, these samples don't exist in expression data, please check sample names"); print(sampleExprInAnnot)}
sampleAnnotInExpr<-cn(ordered_data)[!cn(ordered_data)%in%rn(sampleAnnot)]; if(len(sampleAnnotInExpr)>0){warning("these samples don't exist in sample table, they will be removed from expression data"); print(sampleAnnotInExpr)}
print("Sample names match is ok")

# Create clusters of heatmaps
#exprDE.scaled<-rowScale(ordered_data,center = TRUE,scaled = TRUE)
#quantile.expr<-quantile(unlist(exprDE.scaled),seq(0,1,.01))
#colHA<-colorRamp2(c(quantile.expr[2],0,quantile.expr[100]),c("blue","white","red"))
Filteredmat<-list()
FilteredhclustGeneDE<-list()
FilteredhclustSampleDE<-list()
FilteredHt<-list()
FilteredcolHA <- list()
Filteredha <- list()

mat<-list()
hclustGeneDE<-list()
hclustSampleDE<-list()
Ht<-list()
colHA <- list()
ha <- list()

PVal<- 0.05
FoldChange<-0.58
SampleTokeep<-lire("data/SamplePerCond.tsv")
max<-0
for (i in 1:length(filenames)){
  print(paste("Running heatmap for comparison :",i,"/",length(filenames)))
  file<-lire(filenames[i])
  FilteredData<-ordered_data[row.names(file)[file$Qvalue<PVal&file$AbsoluteAVGLog2Ratio>FoldChange],which(colnames(ordered_data)%in%unlist(SampleTokeep[filenames[i],]))]
  Filteredmat[[i]]<-rowScale(FilteredData,center = TRUE,scaled = TRUE)
  FilteredhclustGeneDE[[i]]<-unsupervisedClustering(Filteredmat[[i]],transpose = F,nboot=30,bootstrap = FALSE,method.dist="pearson")
  FilteredhclustSampleDE[[i]]<-unsupervisedClustering(Filteredmat[[i]],transpose = T,nboot=30,bootstrap = FALSE,method.dist = "euclidean")
  Filteredha[[i]] <- HeatmapAnnotation(Batch = sampleAnnot[FilteredhclustSampleDE[[i]]$labels,"Batch"], Feader = sampleAnnot[FilteredhclustSampleDE[[i]]$labels,"Feader"],Replication = sampleAnnot[FilteredhclustSampleDE[[i]]$labels,"Geminin"],CellType = sampleAnnot[FilteredhclustSampleDE[[i]]$labels,"Tissue"], annotation_legend_param = list(
    Batch = list(
      title = "Batch",
      at = c(1, 2, 3, 4),
      labels = c("1", "2", "3","4"),
      color_bar = "discrete" 
    ),
    Feader = list(
      title = "Feader",
      at = c("DE", "ELB"),
      labels = c("DE", "ELB"),
      color_bar = "discrete" 
    ),
    Replication = list(
      title = "Replication",
      at = c("WithGeminin", "NoGeminin"),
      labels = c("+Gem", "-Gem"),
      color_bar = "discrete"
    ),
    CellType = list(
        title = "Tissue",
        at = c("Sperm", "Progenitor"),
        labels = c("Sperm", "Progenitor"),
        color_bar = "discrete"
    )), col = list(Batch = c("1" = "green", "2" = "red","3" = "blue","4" = "yellow"), Feader = c("DE" = "pink","ELB" = "brown"), Replication = c("WithGeminin" = "grey","NoGeminin" = "black"),CellType = c("Sperm" = "lightgreen","Progenitor" = "cyan")
    ))
  Filteredquantile.expr<-quantile(unlist(Filteredmat[[i]]),seq(0,1,.01))
  FilteredcolHA[[i]]<-colorRamp2(c(Filteredquantile.expr[2],0,Filteredquantile.expr[100]),c("blue","white","red"))
  FilteredHt[[i]]<-Heatmap(Filteredmat[[i]],row_labels = file[row.names(Filteredmat[[i]]),"Genes"], row_names_gp = autoGparFontSizeMatrix(nrow(Filteredmat[[i]])), top_annotation = Filteredha[[i]],
            cluster_rows = FilteredhclustGeneDE[[i]],col = FilteredcolHA[[i]],
            cluster_columns = FilteredhclustSampleDE[[i]],name=paste0(strsplit(strsplit(filenames[i],split = "/")[[1]][2],split = ".t")[[1]][1]," (",nrow(Filteredmat[[i]])," Proteins)"),column_names_gp = autoGparFontSizeMatrix(ncol(Filteredmat[[i]])) )
  
  if(nrow(Filteredmat[[i]]>max)){max<-nrow(Filteredmat[[i]])}
  Data<-ordered_data[row.names(file)[file$Qvalue<PVal&file$AbsoluteAVGLog2Ratio>FoldChange],]
  mat[[i]]<-rowScale(Data,center = TRUE,scaled = TRUE)
  hclustGeneDE[[i]]<-unsupervisedClustering(mat[[i]],transpose = F,nboot=30,bootstrap = FALSE,method.dist="pearson")
  hclustSampleDE[[i]]<-unsupervisedClustering(mat[[i]],transpose = T,nboot=30,bootstrap = FALSE,method.dist = "euclidean")
  ha[[i]] <- HeatmapAnnotation(Batch = sampleAnnot[hclustSampleDE[[i]]$labels,"Batch"], Feader = sampleAnnot[hclustSampleDE[[i]]$labels,"Feader"],Replication = sampleAnnot[hclustSampleDE[[i]]$labels,"Geminin"],CellType = sampleAnnot[hclustSampleDE[[i]]$labels,"Tissue"], annotation_legend_param = list(
    Batch = list(
      title = "Batch",
      at = c(1, 2, 3, 4),
      labels = c("1", "2", "3","4"),
      color_bar = "discrete" 
    ),
    Feader = list(
      title = "Feader",
      at = c("DE", "ELB"),
      labels = c("DE", "ELB"),
      color_bar = "discrete" 
    ),
    Replication = list(
      title = "Replication",
      at = c("WithGeminin", "NoGeminin"),
      labels = c("+Gem", "-Gem"),
      color_bar = "discrete"
    ),
    CellType = list(
      title = "Tissue",
      at = c("Sperm", "Progenitor"),
      labels = c("Sperm", "Progenitor"),
      color_bar = "discrete"
    )), col = list(Batch = c("1" = "green", "2" = "red","3" = "blue","4" = "yellow"), Feader = c("DE" = "pink","ELB" = "brown"), Replication = c("WithGeminin" = "grey","NoGeminin" = "black"),CellType = c("Sperm" = "lightgreen","Progenitor" = "cyan")
    ))
  quantile.expr<-quantile(unlist(mat[[i]]),seq(0,1,.01))
  colHA[[i]]<-colorRamp2(c(quantile.expr[2],0,quantile.expr[100]),c("blue","white","red"))
  Ht[[i]]<-Heatmap(mat[[i]],row_labels = file[row.names(mat[[i]]),"Genes"], row_names_gp = autoGparFontSizeMatrix(nrow(mat[[i]])), top_annotation = ha[[i]],
                   cluster_rows = hclustGeneDE[[i]],col = colHA[[i]],
                   cluster_columns = hclustSampleDE[[i]],name=paste0(strsplit(strsplit(filenames[i],split = "/")[[1]][2],split = ".t")[[1]][1]," (",nrow(mat[[i]])," Proteins)"),column_names_gp = autoGparFontSizeMatrix(ncol(mat[[i]])) )
  ecrire(mat[[i]],paste0("results/MatrixHeatmap",paste0(strsplit(strsplit(filenames[i],split = "/")[[1]][2],split = ".t")[[1]][1]),"_Pval_",PVal,"_FoldChange_",FoldChange,".tsv"))
}

print("Exporting Heatmap into Heatmap.pdf in Figures folder...")
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


print("Summary of DE Protein analysis...")
Distrib<-c()
for (i in 1:length(filenames)){
  Distrib<-c(Distrib,nrow(Filteredmat[[i]]))
}
names(Distrib)<-filenames
pdf(paste0("Figures/Distribution_Barplot_NbrProtDE","_Pval_",PVal,"_FoldChange_",FoldChange,".pdf"),width = 10,height=10)
par(mar=c(15,4,4,2))
barplot(Distrib,main = "Nbr Proteins DE",cex.names = 0.7,las=2,ylim = c(0,8000),names.arg = Distrib)
abline(h = 7314)
dev.off()
ecrire(Distrib,"results/BarplotValues.tsv")

pdf(paste0("Figures/Heatmap_NewMat_AllSamples","_Pval_",PVal,"_FoldChange_",FoldChange,".pdf"),width = 10,height=HeigthScaling)
print(Ht)
dev.off()

pdf(paste0("Figures/Heatmap_NewMat","_Pval_",PVal,"_FoldChange_",FoldChange,".pdf"),width = 10,height=HeigthScaling)
print(FilteredHt)
dev.off()


#### Script for home made Heatmap comparison ####
hclustGeneDE<-list()
Ht<-list()
colHA <- list()
ha <- list()

PVal<- 0.001
FoldChange<-1
SampleTokeep<-lire("data/SamplePerCondHeatmap3.tsv")
DEProt<-c()
for (i in 1:nrow(SampleTokeep)){
  print(paste("Running DE Genes for comparison :",i,"/",row.names(SampleTokeep)[i]))
  file<-lire(row.names(SampleTokeep)[i])
  DEProt<-c(DEProt,row.names(ordered_data[row.names(file)[file$Qvalue<PVal&file$AbsoluteAVGLog2Ratio>FoldChange],]))
}

ProtBackground<-unique(DEProt)
Data<-ordered_data[row.names(file)[which(row.names(ordered_data)%in%ProtBackground)],which(colnames(ordered_data)%in%unlist(unique(SampleTokeep)))]


DataOrder<-Data[,c(13,14,15,16,9,10,12,11,1,4,3,2,5,7,6,8)] #Heatmap 1
DataOrder<-Data[,c(13,14,15,16,9,10,12,11,1,2,5,6,3,4,7,8)] #Heatmap 2 
DataOrder<-Data[,c(1,7,5,4,9,13,11,15,2,8,6,3,10,14,12,16)] #Heatmap 3 
DataOrder<-Data[,c(2,3,9,10,5,7,13,15,1,4,11,12,6,8,14,16)] #Heatmap 4 

mat<-rowScale(DataOrder,center = TRUE,scaled = TRUE)
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


pdf("Figures/HeatmapReplicationELB.pdf",width = 10,height=HeigthScaling)
print(Ht)
dev.off()
ecrire(mat,paste0("results/MatrixHeatmapReplicationELB","_Pval_",PVal,"_FoldChange_",FoldChange,".tsv"))


#### Script for candidates Genes list ####
DEProt<-read.table("ListOfDemandes/list_orthologs_EPI_and_TF_SEGM_DE/list_orthologs_EPI_and_TF_SEGM_DE.txt")
ProtBackground<-unique(DEProt$V1)
Data<-ordered_data[which(row.names(ordered_data)%in%ProtBackground),c(37,38,39,40,33,34,36,35,1,11,5,4,13,19,15,21,2,12,6,3,14,20,16,22)]


c(37,38,39,40,33,34,36,35,1,11,5,4,13,19,15,21) # Sperm prog spermgem- proggem-
c(37,38,39,40,33,34,36,35,1,11,5,4,13,19,15,21,2,12,6,3,14,20,16,22) #Sperm prog spermgem- proggem- spermgem+ proggem+

c(1,11,5,4,13,19,15,21,2,12,6,3,14,20,16,22) # spermgem- proggem- spermgem+ proggem+
c(37,38,39,40,1,11,5,4) #Heatmap sperm vs sperm DE gem-
c(37,38,39,40,1,11,5,4,2,12,6,3)# sperm vs sperm DE gem- vs sperm DE gem+
c(33,34,36,35,13,19,15,21)# prog vs prog DE gem-
c(33,34,36,35,13,19,15,21,14,20,16,22)# prog vs prog DE gem- vs prog DE gem+
c(37,38,39,40,33,34,36,35) # sperm vs prog
i=42
file<-lire(filenames[i])
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

pdf("ListOfDemandes/Heatmap-list_orthologs_EPI_and_TF_SEGM_DE_mail2.pdf",width = 10,height=HeigthScaling)
print(Ht)
dev.off()
ecrire(mat,"ListOfDemandes/Matrix-list_orthologs_EPI_and_TF_SEGM_DE_mail2.tsv")


#### PEA analysis #####
library(pathfindR)
library(filesstrings)
Orthologs<-lire("data/Orthologs.tsv")

Pval<-0.05
database<-"GO-BP"

# Import files to use for the analyze
filenames <- list.files(path="Comparaison", pattern="*.tsv", full.names=TRUE, all.files = TRUE, recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE)
print("Files used for the analysis :")
print(filenames)

i=28
for (i in 28:length(filenames)){
  print(paste("Running heatmap for comparison :",i,"/",length(filenames)))
  file<-lire(filenames[i])
  Input<-file[,c(8,2,5)]  #Human gene name, logFC, PVal
  Input<-Input[which(row.names(Input)%in%row.names(Orthologs)),]
  Input_ordered<-Input[row.names(Orthologs),]
  Input_ordered$Genes<-Orthologs$HumanName
  output_all<-run_pathfindR(Input_ordered,p_val_threshold=Pval,gene_sets=database)
  output_pos<-run_pathfindR(Input_ordered[which(Input_ordered$AVGLog2Ratio>0),],p_val_threshold=Pval,gene_sets=database)
  output_neg<-run_pathfindR(Input_ordered[which(Input_ordered$AVGLog2Ratio<0),],p_val_threshold=Pval,gene_sets=database)
  pdf(paste0("results/GoEnrichment",paste0(strsplit(strsplit(filenames[i],split = "/")[[1]][2],split = ".t")[[1]][1]),"_Pval_",Pval,"_database_",database,".pdf"),width = 10,height=10)
  EnrichmentGraph<-enrichment_chart(
    result_df = output_all,
    top_terms = 30
  )
  
  EnrichmentGraphPos<-enrichment_chart(
    result_df = output_pos,
    top_terms = 30
  )
  
  EnrichmentGraphNeg<-enrichment_chart(
    result_df = output_neg,
    top_terms = 30
  )
  
  #input_processed <- input_processing(Input_ordered)
  
  #BioGrid_Graph<-visualize_terms(result_df = output_all, input_processed = input_processed, hsa_KEGG = FALSE, pin_name_path = "Biogrid")
  
  #dir.create(paste0("term_visualizations/",paste0(strsplit(strsplit(filenames[i],split = "/")[[1]][2],split = ".t")[[1]][1]),"_Pval_",PVal,"_database_",database))
  #file.move(paste0("term_visualizations/",list.files(path="term_visualizations", pattern=".png")), paste0("term_visualizations/",paste0(strsplit(strsplit(filenames[i],split = "/")[[1]][2],split = ".t")[[1]][1]),"_Pval_",PVal,"_database_",database))
  
  HeatmapGraph<-term_gene_heatmap(result_df = output_all, genes_df = Input_ordered)
  NonOrientedGraph<-term_gene_graph(result_df = output_all, use_description = TRUE)
  UpsetPlot<-UpSet_plot(result_df = output_all, genes_df = Input_ordered)
  print(EnrichmentGraph)
  print(EnrichmentGraphPos)
  print(EnrichmentGraphNeg)
  print(HeatmapGraph)
  print(NonOrientedGraph)
  print(UpsetPlot)
  dev.off()
}


#### Transcriptomic combinatory ####
row.names(ordered_data)<-file[row.names(ordered_data),"Genes"]
PrProtValues<-as.matrix(apply(ordered_data[,c(33,34,36,35)],1,mean))
SpProtValues<-as.matrix(apply(ordered_data[,c(37,38,39,40)],1,mean))
colnames(PrProtValues)<-c("ProgProt")
colnames(SpProtValues)<-c("SpProt")

Temp<-exprW[which(row.names(exprW)%in%file$Genes),]
file<-file[order(file$Genes),]
row.names(Temp)<-file$ProteinGroups[which(file$Genes%in%row.names(Temp))]
RNAValues<-as.matrix(log2(apply(Temp[,c(1,2,3)],1,mean)/apply(Temp[,c(4,5)],1,mean)))
RNAValues<-RNAValues[order(row.names(RNAValues)),]

PrProtValues<-PrProtValues[which(row.names(PrProtValues)%in%row.names(FoldchangeMatrix)),]
SpProtValues<-SpProtValues[which(row.names(SpProtValues)%in%row.names(FoldchangeMatrix)),]
FoldchangeMatrix<-cbind(FoldchangeMatrix,PrProtValues,SpProtValues)

setwd("~/These/Analyses/NolwennProject")

i=37
file<-lire(filenames[i])
setwd("~/These/Analyses/NolwennProject/RNA-Seq/")
FoldchangeMatrix<-lire("results/MatrixVolcanoPlotPerso.tsv")
logFCthreshold<-0.58 #Absolute Log2(Fold-Change) threshold (if logFCthreshold=1, gene is differentially expressed if expressed 2 time more or less between folds)
ProtDE<-row.names(lire(x = "../ListOfDemandes/list_ortholog_up_down_sperm_vs_prog/list_ortholog_up_down_sperm_vs_prog.txt"))

pdf(paste0("figs/VolcanoPlotColor.pdf"),width=10,height=10)
plot(FoldchangeMatrix[,2],FoldchangeMatrix[,1],type="n",ylab="Progenitor RNA-Seq",xlab="log2FC(Sp/Pr) Protein",col="black",main=paste0("Volcano plot"))
abline(v = -logFCthreshold)
abline(v = logFCthreshold)
text(FoldchangeMatrix[,2],FoldchangeMatrix[,1],labels = file[row.names(FoldchangeMatrix),"Genes"],cex=.2,col="black")
text(FoldchangeMatrix[which(row.names(FoldchangeMatrix)%in%ProtDE),2],FoldchangeMatrix[which(row.names(FoldchangeMatrix)%in%ProtDE),1],labels = file[row.names(FoldchangeMatrix[which(row.names(FoldchangeMatrix)%in%ProtDE),]),"Genes"],cex=.2,col="red")
dev.off()

pdf(paste0("figs/VolcanoPlotAutoProtValuesColorAutoFC.pdf"),width=10,height=10)
plot(FoldchangeMatrix[,2],FoldchangeMatrix[,1],type="n",ylab="Progenitor RNA-Seq",xlab="AVGLog2Ratio",col="black",main=paste0("Volcano plot"))
abline(v = -logFCthreshold)
abline(v = logFCthreshold)
text(file[row.names(FoldchangeMatrix),"AVGLog2Ratio"],FoldchangeMatrix[,1],labels = file[row.names(FoldchangeMatrix),"Genes"],cex=.2,col="black")
text(file[row.names(FoldchangeMatrix[which(row.names(FoldchangeMatrix)%in%ProtDE),]),"AVGLog2Ratio"],FoldchangeMatrix[which(row.names(FoldchangeMatrix)%in%ProtDE),1],labels = file[row.names(FoldchangeMatrix[which(row.names(FoldchangeMatrix)%in%ProtDE),]),"Genes"],cex=.2,col="red")
dev.off()

pdf(paste0("figs/VolcanoPlotRNAvsProtValues.pdf"),width=10,height=10)
plot(FoldchangeMatrix[,3],FoldchangeMatrix[,1],type="n",ylab="Progenitor RNA-Seq",xlab="Prog Protein",col="black",main=paste0("Volcano plot"))
text(FoldchangeMatrix[,3],FoldchangeMatrix[,1],labels = file[row.names(FoldchangeMatrix),"Genes"],cex=.2,col="black")
text(FoldchangeMatrix[which(row.names(FoldchangeMatrix)%in%ProtDE),3],FoldchangeMatrix[which(row.names(FoldchangeMatrix)%in%ProtDE),1],labels = file[row.names(FoldchangeMatrix[which(row.names(FoldchangeMatrix)%in%ProtDE),]),"Genes"],cex=.2,col="red")

plot(FoldchangeMatrix[,4],FoldchangeMatrix[,1],type="n",ylab="Progenitor RNA-Seq",xlab="Sperm Protein",col="black",main=paste0("Volcano plot"))
text(FoldchangeMatrix[,4],FoldchangeMatrix[,1],labels = file[row.names(FoldchangeMatrix),"Genes"],cex=.2,col="black")
text(FoldchangeMatrix[which(row.names(FoldchangeMatrix)%in%ProtDE),4],FoldchangeMatrix[which(row.names(FoldchangeMatrix)%in%ProtDE),1],labels = file[row.names(FoldchangeMatrix[which(row.names(FoldchangeMatrix)%in%ProtDE),]),"Genes"],cex=.2,col="red")

dev.off()

ecrire(x=FoldchangeMatrix,file="results/MatrixVolcanoPlotPerso.tsv")

pdf(paste0("figs/FACSlikeLessStringeant.pdf"),width=10,height=10)
ToPlot<-4
hist(FoldchangeMatrix[,3]) # 6 et 13
hist(FoldchangeMatrix[,1]) # 9 et 16
DoubleNeg<-FoldchangeMatrix[which(FoldchangeMatrix$RNAValues<10&FoldchangeMatrix$PrProtValues<7),]
plot(DoubleNeg[,ToPlot],DoubleNeg[,1],type="n",ylab="Progenitor RNA-Seq",xlab="Sperm Protein",col="black",main=paste0("Volcano plot double neg"),xlim = c(0,20),ylim=c(0,20))
text(DoubleNeg[,ToPlot],DoubleNeg[,1],labels = file[row.names(DoubleNeg),"Genes"],cex=.2,col="black")
text(DoubleNeg[which(row.names(DoubleNeg)%in%ProtDE),ToPlot],DoubleNeg[which(row.names(DoubleNeg)%in%ProtDE),1],labels = file[row.names(DoubleNeg[which(row.names(DoubleNeg)%in%ProtDE),]),"Genes"],cex=.2,col="red")

DoublePos<-FoldchangeMatrix[which(FoldchangeMatrix$RNAValues>14&FoldchangeMatrix$PrProtValues>12),]
plot(DoublePos[,ToPlot],DoublePos[,1],type="n",ylab="Progenitor RNA-Seq",xlab="Sperm Protein",col="black",main=paste0("Volcano plot double pos"),xlim = c(0,20),ylim=c(0,20))
text(DoublePos[,ToPlot],DoublePos[,1],labels = file[row.names(DoublePos),"Genes"],cex=.2,col="black")
text(DoublePos[which(row.names(DoublePos)%in%ProtDE),ToPlot],DoublePos[which(row.names(DoublePos)%in%ProtDE),1],labels = file[row.names(DoublePos[which(row.names(DoublePos)%in%ProtDE),]),"Genes"],cex=.2,col="red")

HautGauche<-FoldchangeMatrix[which(FoldchangeMatrix$RNAValues>14&FoldchangeMatrix$PrProtValues<7),]
plot(HautGauche[,ToPlot],HautGauche[,1],type="n",ylab="Progenitor RNA-Seq",xlab="Sperm Protein",col="black",main=paste0("Volcano plot haut gauche"),xlim = c(0,20),ylim=c(0,20))
text(HautGauche[,ToPlot],HautGauche[,1],labels = file[row.names(HautGauche),"Genes"],cex=.2,col="black")
text(HautGauche[which(row.names(HautGauche)%in%ProtDE),ToPlot],HautGauche[which(row.names(HautGauche)%in%ProtDE),1],labels = file[row.names(HautGauche[which(row.names(HautGauche)%in%ProtDE),]),"Genes"],cex=.2,col="red")

BasDroite<-FoldchangeMatrix[which(FoldchangeMatrix$RNAValues<10&FoldchangeMatrix$PrProtValues>12),]
plot(BasDroite[,ToPlot],BasDroite[,1],type="n",ylab="Progenitor RNA-Seq",xlab="Sperm Protein",col="black",main=paste0("Volcano plot bas droite"),xlim = c(0,20),ylim=c(0,20))
text(BasDroite[,ToPlot],BasDroite[,1],labels = file[row.names(BasDroite),"Genes"],cex=.2,col="black")
text(BasDroite[which(row.names(BasDroite)%in%ProtDE),ToPlot],BasDroite[which(row.names(BasDroite)%in%ProtDE),1],labels = file[row.names(BasDroite[which(row.names(BasDroite)%in%ProtDE),]),"Genes"],cex=.2,col="red")

dev.off()

pdf(paste0("figs/FACSlikeLessStringeantFC.pdf"),width=10,height=10)
ToPlot<-2
hist(FoldchangeMatrix[,3]) # 6 et 13
hist(FoldchangeMatrix[,1]) # 9 et 16
DoubleNeg<-FoldchangeMatrix[which(FoldchangeMatrix$RNAValues<10&FoldchangeMatrix$PrProtValues<7),]
plot(DoubleNeg[,ToPlot],DoubleNeg[,1],type="n",ylab="Progenitor RNA-Seq",xlab="log2FC(Sp/Pr) Protein",col="black",main=paste0("Volcano plot double neg"),xlim = c(-2,2),ylim=c(0,20))
text(DoubleNeg[,ToPlot],DoubleNeg[,1],labels = file[row.names(DoubleNeg),"Genes"],cex=.2,col="black")
text(DoubleNeg[which(row.names(DoubleNeg)%in%ProtDE),ToPlot],DoubleNeg[which(row.names(DoubleNeg)%in%ProtDE),1],labels = file[row.names(DoubleNeg[which(row.names(DoubleNeg)%in%ProtDE),]),"Genes"],cex=.2,col="red")

DoublePos<-FoldchangeMatrix[which(FoldchangeMatrix$RNAValues>14&FoldchangeMatrix$PrProtValues>12),]
plot(DoublePos[,ToPlot],DoublePos[,1],type="n",ylab="Progenitor RNA-Seq",xlab="log2FC(Sp/Pr) Protein",col="black",main=paste0("Volcano plot double pos"),xlim = c(-2,2),ylim=c(0,20))
text(DoublePos[,ToPlot],DoublePos[,1],labels = file[row.names(DoublePos),"Genes"],cex=.2,col="black")
text(DoublePos[which(row.names(DoublePos)%in%ProtDE),ToPlot],DoublePos[which(row.names(DoublePos)%in%ProtDE),1],labels = file[row.names(DoublePos[which(row.names(DoublePos)%in%ProtDE),]),"Genes"],cex=.2,col="red")

HautGauche<-FoldchangeMatrix[which(FoldchangeMatrix$RNAValues>14&FoldchangeMatrix$PrProtValues<7),]
plot(HautGauche[,ToPlot],HautGauche[,1],type="n",ylab="Progenitor RNA-Seq",xlab="log2FC(Sp/Pr) Protein",col="black",main=paste0("Volcano plot haut gauche"),xlim = c(-2,2),ylim=c(0,20))
text(HautGauche[,ToPlot],HautGauche[,1],labels = file[row.names(HautGauche),"Genes"],cex=.2,col="black")
text(HautGauche[which(row.names(HautGauche)%in%ProtDE),ToPlot],HautGauche[which(row.names(HautGauche)%in%ProtDE),1],labels = file[row.names(HautGauche[which(row.names(HautGauche)%in%ProtDE),]),"Genes"],cex=.2,col="red")

BasDroite<-FoldchangeMatrix[which(FoldchangeMatrix$RNAValues<10&FoldchangeMatrix$PrProtValues>12),]
plot(BasDroite[,ToPlot],BasDroite[,1],type="n",ylab="Progenitor RNA-Seq",xlab="log2FC(Sp/Pr) Protein",col="black",main=paste0("Volcano plot bas droite"),xlim = c(-2,2),ylim=c(0,20))
text(BasDroite[,ToPlot],BasDroite[,1],labels = file[row.names(BasDroite),"Genes"],cex=.2,col="black")
text(BasDroite[which(row.names(BasDroite)%in%ProtDE),ToPlot],BasDroite[which(row.names(BasDroite)%in%ProtDE),1],labels = file[row.names(BasDroite[which(row.names(BasDroite)%in%ProtDE),]),"Genes"],cex=.2,col="red")

dev.off()

pdf(paste0("figs/FACSlikeLessStringeantAutoFC.pdf"),width=10,height=10)
ToPlot<-2
hist(FoldchangeMatrix[,3]) # 6 et 13
hist(FoldchangeMatrix[,1]) # 9 et 16
DoubleNeg<-FoldchangeMatrix[which(FoldchangeMatrix$RNAValues<10&FoldchangeMatrix$PrProtValues<7),]
plot(DoubleNeg[,ToPlot],DoubleNeg[,1],type="n",ylab="Progenitor RNA-Seq",xlab="log2FC(Sp/Pr) Protein",col="black",main=paste0("Volcano plot double neg"),xlim = c(-2,2),ylim=c(0,20))
text(file[row.names(DoubleNeg),"AVGLog2Ratio"],DoubleNeg[,1],labels = file[row.names(DoubleNeg),"Genes"],cex=.2,col="black")
text(file[row.names(DoubleNeg[which(row.names(DoubleNeg)%in%ProtDE),]),"AVGLog2Ratio"],DoubleNeg[which(row.names(DoubleNeg)%in%ProtDE),1],labels = file[row.names(DoubleNeg[which(row.names(DoubleNeg)%in%ProtDE),]),"Genes"],cex=.2,col="red")

DoublePos<-FoldchangeMatrix[which(FoldchangeMatrix$RNAValues>14&FoldchangeMatrix$PrProtValues>12),]
plot(DoublePos[,ToPlot],DoublePos[,1],type="n",ylab="Progenitor RNA-Seq",xlab="log2FC(Sp/Pr) Protein",col="black",main=paste0("Volcano plot double pos"),xlim = c(-2,2),ylim=c(0,20))
text(file[row.names(DoublePos),"AVGLog2Ratio"],DoublePos[,1],labels = file[row.names(DoublePos),"Genes"],cex=.2,col="black")
text(file[row.names(DoublePos[which(row.names(DoublePos)%in%ProtDE),]),"AVGLog2Ratio"],DoublePos[which(row.names(DoublePos)%in%ProtDE),1],labels = file[row.names(DoublePos[which(row.names(DoublePos)%in%ProtDE),]),"Genes"],cex=.2,col="red")

HautGauche<-FoldchangeMatrix[which(FoldchangeMatrix$RNAValues>14&FoldchangeMatrix$PrProtValues<7),]
plot(HautGauche[,ToPlot],HautGauche[,1],type="n",ylab="Progenitor RNA-Seq",xlab="log2FC(Sp/Pr) Protein",col="black",main=paste0("Volcano plot haut gauche"),xlim = c(-2,2),ylim=c(0,20))
text(file[row.names(HautGauche),"AVGLog2Ratio"],HautGauche[,1],labels = file[row.names(HautGauche),"Genes"],cex=.2,col="black")
text(file[row.names(HautGauche[which(row.names(HautGauche)%in%ProtDE),]),"AVGLog2Ratio"],HautGauche[which(row.names(HautGauche)%in%ProtDE),1],labels = file[row.names(HautGauche[which(row.names(HautGauche)%in%ProtDE),]),"Genes"],cex=.2,col="red")

BasDroite<-FoldchangeMatrix[which(FoldchangeMatrix$RNAValues<10&FoldchangeMatrix$PrProtValues>12),]
plot(BasDroite[,ToPlot],BasDroite[,1],type="n",ylab="Progenitor RNA-Seq",xlab="log2FC(Sp/Pr) Protein",col="black",main=paste0("Volcano plot bas droite"),xlim = c(-2,2),ylim=c(0,20))
text(file[row.names(BasDroite),"AVGLog2Ratio"],BasDroite[,1],labels = file[row.names(BasDroite),"Genes"],cex=.2,col="black")
text(file[row.names(BasDroite[which(row.names(BasDroite)%in%ProtDE),]),"AVGLog2Ratio"],BasDroite[which(row.names(BasDroite)%in%ProtDE),1],labels = file[row.names(BasDroite[which(row.names(BasDroite)%in%ProtDE),]),"Genes"],cex=.2,col="red")

dev.off()
