library(dplyr)
library(Seurat)
library(patchwork)

BCP002_Total_3GEX<- as.data.frame(Read10X(data.dir = "c:/Users/xjmik/Downloads/E-MTAB-9005.processed.2/BCP002_Total_3GEX_filtered_feature_bc_matrix/"))
BCP003_Total<- as.data.frame(Read10X(data.dir = "c:/Users/xjmik/Downloads/E-MTAB-9005.processed.2/BCP003_Total_5GEX_filtered_feature_bc_matrix/"))
BCP003_MBC<-as.data.frame(Read10X(data.dir = "c:/Users/xjmik/Downloads/E-MTAB-9005.processed.2/BCP003_MBC_5GEX_filtered_feature_bc_matrix/"))
BCP004_Total<-as.data.frame(Read10X(data.dir = "c:/Users/xjmik/Downloads/E-MTAB-9005.processed.2/BCP004_Total_5GEX_filtered_feature_bc_matrix/"))
BCP004_MBC<-as.data.frame(Read10X(data.dir = "c:/Users/xjmik/Downloads/E-MTAB-9005.processed.2/BCP004_MBC_5GEX_filtered_feature_bc_matrix/"))
BCP005_Total<-as.data.frame(Read10X(data.dir = "c:/Users/xjmik/Downloads/E-MTAB-9005.processed.2/BCP005_Total_5GEX_filtered_feature_bc_matrix/"))
BCP005_MBC<-as.data.frame(Read10X(data.dir = "c:/Users/xjmik/Downloads/E-MTAB-9005.processed.2/BCP005_IgMnegMBC_5GEX_filtered_feature_bc_matrix/"))
BCP006_Total<-as.data.frame(Read10X(data.dir = "c:/Users/xjmik/Downloads/E-MTAB-9005.processed.2/BCP006_Total_5GEX_filtered_feature_bc_matrix/"))
BCP006_MBC<-as.data.frame(Read10X(data.dir = "c:/Users/xjmik/Downloads/E-MTAB-9005.processed.2/BCP006_IgMnegMBC_5GEX_filtered_feature_bc_matrix/"))
BCP008_Total<-as.data.frame(Read10X(data.dir = "c:/Users/xjmik/Downloads/E-MTAB-9005.processed.2/BCP008_Total_5GEX_filtered_feature_bc_matrix/"))
BCP008_MBC<-as.data.frame(Read10X(data.dir = "c:/Users/xjmik/Downloads/E-MTAB-9005.processed.2/BCP008_IgMneg_MBC_5GEX_filtered_feature_bc_matrix/"))
BCP009_Total<-as.data.frame(Read10X(data.dir = "c:/Users/xjmik/Downloads/E-MTAB-9005.processed.2/BCP009_Total_5GEX_filtered_feature_bc_matrix/"))
BCP009_MBC<-as.data.frame(Read10X(data.dir = "c:/Users/xjmik/Downloads/E-MTAB-9005.processed.2/BCP009_IgMneg_MBC_5GEX_filtered_feature_bc_matrix/"))

for (i in 1:length(colnames(BCP002_Total_3GEX))) {
  colnames(BCP002_Total_3GEX)[i] <- paste(colnames(BCP002_Total_3GEX)[i],"BCP2_Total",i,sep = "-")
}

for (i in 1:length(colnames(BCP003_MBC))) {
  colnames(BCP003_MBC)[i] <- paste(colnames(BCP003_MBC)[i],"BCP3_MBC",i,sep = "-")
}

for (i in 1:length(colnames(BCP003_Total))) {
  colnames(BCP003_Total)[i] <- paste(colnames(BCP003_Total)[i],"BCP3_Total",i,sep = "-")
}

for (i in 1:length(colnames(BCP004_MBC))) {
  colnames(BCP004_MBC)[i] <- paste(colnames(BCP004_MBC)[i],"BCP4_MBC",i,sep = "-")
}

for (i in 1:length(colnames(BCP004_Total))) {
  colnames(BCP004_Total)[i] <- paste(colnames(BCP004_Total)[i],"BCP4_Total",i,sep = "-")
}

for (i in 1:length(colnames(BCP005_MBC))) {
  colnames(BCP005_MBC)[i] <- paste(colnames(BCP005_MBC)[i],"BCP5_MBC",i,sep = "-")
}

for (i in 1:length(colnames(BCP005_Total))) {
  colnames(BCP005_Total)[i] <- paste(colnames(BCP005_Total)[i],"BCP5_Total",i,sep = "-")
}

for (i in 1:length(colnames(BCP006_MBC))) {
  colnames(BCP006_MBC)[i] <- paste(colnames(BCP006_MBC)[i],"BCP6_MBC",i,sep = "-")
}

for (i in 1:length(colnames(BCP006_Total))) {
  colnames(BCP006_Total)[i] <- paste(colnames(BCP006_Total)[i],"BCP6_Total",i,sep = "-")
}

for (i in 1:length(colnames(BCP008_MBC))) {
  colnames(BCP008_MBC)[i] <- paste(colnames(BCP008_MBC)[i],"BCP8_MBC",i,sep = "-")
}

for (i in 1:length(colnames(BCP008_Total))) {
  colnames(BCP008_Total)[i] <- paste(colnames(BCP008_Total)[i],"BCP8_Total",i,sep = "-")
}

for (i in 1:length(colnames(BCP009_MBC))) {
  colnames(BCP009_MBC)[i] <- paste(colnames(BCP009_MBC)[i],"BCP9_MBC",i,sep = "-")
}

for (i in 1:length(colnames(BCP009_Total))) {
  colnames(BCP009_Total)[i] <- paste(colnames(BCP009_Total)[i],"BCP9_Total",i,sep = "-")
}
BCP002_Total_3GEX.metadata<-data.frame(colnames(BCP002_Total_3GEX),rep("BCP2",length(colnames(BCP002_Total_3GEX))),rep("Total",length(colnames(BCP002_Total_3GEX))))
BCP003_MBC.metadata<-data.frame(colnames(BCP003_MBC),rep("BCP3",length(colnames(BCP003_MBC))),rep("MBC",length(colnames(BCP003_MBC))))
BCP003_Total.metadata<-data.frame(colnames(BCP003_Total),rep("BCP3",length(colnames(BCP003_Total))),rep("Total",length(colnames(BCP003_Total))))
BCP004_MBC.metadata<-data.frame(colnames(BCP004_MBC),rep("BCP4",length(colnames(BCP004_MBC))),rep("MBC",length(colnames(BCP004_MBC))))
BCP004_Total.metadata<-data.frame(colnames(BCP004_Total),rep("BCP4",length(colnames(BCP004_Total))),rep("Total",length(colnames(BCP004_Total))))
BCP005_MBC.metadata<-data.frame(colnames(BCP005_MBC),rep("BCP5",length(colnames(BCP005_MBC))),rep("MBC",length(colnames(BCP005_MBC))))
BCP005_Total.metadata<-data.frame(colnames(BCP005_Total),rep("BCP5",length(colnames(BCP005_Total))),rep("Total",length(colnames(BCP005_Total))))
BCP006_MBC.metadata<-data.frame(colnames(BCP006_MBC),rep("BCP6",length(colnames(BCP006_MBC))),rep("MBC",length(colnames(BCP006_MBC))))
BCP006_Total.metadata<-data.frame(colnames(BCP006_Total),rep("BCP6",length(colnames(BCP006_Total))),rep("Total",length(colnames(BCP006_Total))))
BCP008_MBC.metadata<-data.frame(colnames(BCP008_MBC),rep("BCP8",length(colnames(BCP008_MBC))),rep("MBC",length(colnames(BCP008_MBC))))
BCP008_Total.metadata<-data.frame(colnames(BCP008_Total),rep("BCP8",length(colnames(BCP008_Total))),rep("Total",length(colnames(BCP008_Total))))
BCP009_MBC.metadata<-data.frame(colnames(BCP009_MBC),rep("BCP9",length(colnames(BCP009_MBC))),rep("MBC",length(colnames(BCP009_MBC))))
BCP009_Total.metadata<-data.frame(colnames(BCP009_Total),rep("BCP9",length(colnames(BCP009_Total))),rep("Total",length(colnames(BCP009_Total))))

colnames(BCP002_Total_3GEX.metadata)<-c("barcode","patient","group")
colnames(BCP003_MBC.metadata)<-c("barcode","patient","group")
colnames(BCP003_Total.metadata)<-c("barcode","patient","group")
colnames(BCP004_MBC.metadata)<-c("barcode","patient","group")
colnames(BCP004_Total.metadata)<-c("barcode","patient","group")
colnames(BCP005_MBC.metadata)<-c("barcode","patient","group")
colnames(BCP005_Total.metadata)<-c("barcode","patient","group")
colnames(BCP006_MBC.metadata)<-c("barcode","patient","group")
colnames(BCP006_Total.metadata)<-c("barcode","patient","group")
colnames(BCP008_MBC.metadata)<-c("barcode","patient","group")
colnames(BCP008_Total.metadata)<-c("barcode","patient","group")
colnames(BCP009_MBC.metadata)<-c("barcode","patient","group")
colnames(BCP009_Total.metadata)<-c("barcode","patient","group")

pbmc.metadata<-rbind(BCP002_Total_3GEX.metadata,BCP003_MBC.metadata,BCP003_Total.metadata,BCP004_MBC.metadata,BCP004_Total.metadata,BCP005_MBC.metadata,BCP005_Total.metadata,BCP006_MBC.metadata,BCP006_Total.metadata,BCP008_MBC.metadata,BCP008_Total.metadata,BCP009_MBC.metadata,BCP009_Total.metadata)
remove(BCP002_Total_3GEX.metadata,BCP003_MBC.metadata,BCP003_Total.metadata,BCP004_MBC.metadata,BCP004_Total.metadata,BCP005_MBC.metadata,BCP005_Total.metadata,BCP006_MBC.metadata,BCP006_Total.metadata,BCP008_MBC.metadata,BCP008_Total.metadata,BCP009_MBC.metadata,BCP009_Total.metadata)
pbmc.data<-cbind(BCP002_Total_3GEX,BCP003_MBC,BCP003_Total,BCP004_MBC,BCP004_Total,BCP005_MBC,BCP005_Total,BCP006_MBC,BCP006_Total,BCP008_MBC,BCP008_Total,BCP009_MBC,BCP009_Total)
remove(BCP002_Total_3GEX,BCP003_MBC,BCP003_Total,BCP004_MBC,BCP004_Total,BCP005_MBC,BCP005_Total,BCP006_MBC,BCP006_Total,BCP008_MBC,BCP008_Total,BCP009_MBC,BCP009_Total)
remove(i)

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200,meta.data = pbmc.metadata )
#remove(pbmc.data,pbmc.metadata)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
counts<-quantile(pbmc@meta.data$nCount_RNA,c(0.025,0.975))
Features<-quantile(pbmc@meta.data$nFeature_RNA,c(0.025,0.975))
percentmt<-quantile(pbmc@meta.data$percent.mt,c(0.025,0.975))
pbmc <- subset(pbmc, subset = nFeature_RNA > 397 & nFeature_RNA < 3895 & percent.mt < 26.112 & nCount_RNA < 18980 & nCount_RNA >701 & percent.mt > 0.175)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc1 <- ScaleData(pbmc,features = all.genes)
pbmc2 <- ScaleData(pbmc,features = all.genes,model.use = 'poison')
# samle as pbmc2 pbmc3 <- ScaleData(pbmc,features = all.genes,model.use = 'negbinom')
pbmc3 <- ScaleData(pbmc,features = all.genes,vars.to.regress = "percent.mt")
pbmc4 <- ScaleData(pbmc,features = all.genes,vars.to.regress = "nFeature_RNA")
pbmc5 <- ScaleData(pbmc,features = all.genes,vars.to.regress = "nCount_RNA")
pbmc6 <- ScaleData(pbmc,features = all.genes,vars.to.regress = c("percent.mt","nCount_RNA","nFeature_RNA"))

