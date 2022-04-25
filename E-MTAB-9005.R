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
metadata<-read.table("c:/Users/xjmik/Downloads/E-MTAB-9005.processed.1/CellTypeMetaData.txt",sep = "\t",header = TRUE,row.names = 1)
a<-as.data.frame(rownames(metadata))
colnames(a)<-"metadata"
library(stringr)
pbmc.metadata[,c("a","b","c","d")]<-str_split_fixed(pbmc.metadata$barcode,"-",4)
a[,c("e","f","g")]<-str_split_fixed(a$metadata,"_",3)

rownames(pbmc.metadata)<-pbmc.metadata$barcode
rownames(a)<-rownames(metadata)
metadata<-cbind(a,metadata)
remove(a)

pbmc.metadata$h<-paste0(pbmc.metadata$c,"_",pbmc.metadata$a)
rownames(pbmc.metadata)<-pbmc.metadata$h
pbmc.metadata1<-pbmc.metadata[rownames(metadata),]
pbmc.metadata1<-cbind(pbmc.metadata1,metadata)
rownames(pbmc.metadata1)<-pbmc.metadata1$barcode
pbmc.data1<-pbmc.data[,rownames(pbmc.metadata1)]
remove(metadata)
rownames(pbmc.metadata)<-pbmc.metadata$barcode

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200,meta.data = pbmc.metadata )
pbmc1 <- CreateSeuratObject(counts = pbmc.data1, project = "pbmc3k", min.cells = 3, min.features = 200,meta.data = pbmc.metadata1 )
remove(pbmc.data,pbmc.data1,pbmc.metadata,pbmc.metadata1)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc1[["percent.mt"]] <- PercentageFeatureSet(pbmc1, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
VlnPlot(pbmc1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
pbmc1 <- subset(pbmc1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)
pbmc <- NormalizeData(pbmc)
pbmc1 <- NormalizeData(pbmc1)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc1 <- FindVariableFeatures(pbmc1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
all.genes <- rownames(pbmc1)
pbmc1 <- ScaleData(pbmc1, features = all.genes)
remove(all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc1 <- RunPCA(pbmc1, features = VariableFeatures(object = pbmc1))
pbmc_1 <- RunPCA(pbmc, features = rownames(pbmc))
pbmc1_1 <- RunPCA(pbmc1, features = rownames(pbmc1))

ElbowPlot(pbmc,ndims = 50)
ElbowPlot(pbmc_1,ndims = 50)
ElbowPlot(pbmc1,ndims = 50)
ElbowPlot(pbmc1_1,ndims = 50)

pbmc <- FindNeighbors(pbmc, dims = 1:12)
pbmc <- FindClusters(pbmc, resolution = 0.8)
pbmc_1 <- FindNeighbors(pbmc_1, dims = 1:20)
pbmc_1 <- FindClusters(pbmc_1, resolution = 0.8)
pbmc1 <- FindNeighbors(pbmc1, dims = 1:20)
pbmc1 <- FindClusters(pbmc1, resolution = 0.8)
pbmc1_1 <- FindNeighbors(pbmc1_1, dims = 1:20)
pbmc1_1 <- FindClusters(pbmc1_1, resolution = 0.8)

pbmc <- RunUMAP(pbmc, dims = 1:12)
pbmc <- RunTSNE(pbmc, dims = 1:12)
pbmc_1 <- RunUMAP(pbmc_1, dims = 1:20)
pbmc_1 <- RunTSNE(pbmc_1, dims = 1:20)
pbmc1 <- RunUMAP(pbmc1, dims = 1:20)
pbmc1 <- RunTSNE(pbmc1, dims = 1:20)
pbmc1_1 <- RunUMAP(pbmc1_1, dims = 1:20)
pbmc1_1 <- RunTSNE(pbmc1_1, dims = 1:20)

DimPlot(pbmc, reduction = "umap")
DimPlot(pbmc_1, reduction = "umap")
DimPlot(pbmc1, reduction = "umap")
DimPlot(pbmc1_1, reduction = "umap")

DimPlot(pbmc, reduction = "tsne")
DimPlot(pbmc_1, reduction = "tsne")
DimPlot(pbmc1, reduction = "tsne")
DimPlot(pbmc1_1, reduction = "tsne")

Idents(pbmc1)<-pbmc1@meta.data$Lineage
DimPlot(pbmc1, reduction = "umap")
DimPlot(pbmc1, reduction = "tsne")

Idents(pbmc1_1)<-pbmc1_1@meta.data$Lineage
DimPlot(pbmc1_1, reduction = "umap")
DimPlot(pbmc1, reduction = "tsne")

Idents(pbmc1)<-pbmc1@meta.data$CellType
DimPlot(pbmc1, reduction = "umap",label = TRUE)
DimPlot(pbmc1, reduction = "tsne",label = TRUE)

Idents(pbmc1_1)<-pbmc1_1@meta.data$CellType
DimPlot(pbmc1_1, reduction = "umap",label = TRUE)
DimPlot(pbmc1_1, reduction = "tsne",label = TRUE)

VlnPlot(pbmc1,features = c("MYC","IRF4","KDM6B","LY75"),pt.size = 0,ncol = 2)
VlnPlot(pbmc1_1,features = c("MYC","IRF4","KDM6B","LY75"),pt.size = 0,ncol = 2)

FeaturePlot(pbmc1,features = c("MYC","IRF4","KDM6B","LY75"),label = TRUE)
FeaturePlot(pbmc1,features = c("MYC","IRF4","KDM6B","LY75"),label = TRUE,reduction = "tsne")
FeaturePlot(pbmc1_1,features = c("MYC","IRF4","KDM6B","LY75"),label = TRUE)
FeaturePlot(pbmc1_1,features = c("MYC","IRF4","KDM6B","LY75"),label = TRUE,reduction = "tsne")

library(future)
plan("multiprocess", workers = 20)
pbmc1_Celltype.markers<-FindAllMarkers(pbmc1,only.pos = TRUE,logfc.threshold = 0,min.pct = 0)
pbmc1_1_Celltype.markers<-FindAllMarkers(pbmc1_1,only.pos = TRUE,logfc.threshold = 0,min.pct = 0)
write.table(pbmc1_Celltype.markers,"pbmc1_Celltype.markers.txt",sep = "\t")
write.table(pbmc1_1_Celltype.markers,"pbmc1_1_Celltype.markers.txt",sep = "\t")
remove(pbmc1_1_Celltype.markers,pbmc1_Celltype.markers)

pbmc1_Celltype_averagecluster<-AverageExpression(pbmc1,assays = "RNA")
pbmc1_1_Celltype_averagecluster<-AverageExpression(pbmc1_1,assays = "RNA")
write.table(pbmc1_Celltype_averagecluster$RNA,"pbmc1_Celltype_averagecluster.txt",sep = "\t")
write.table(pbmc1_1_Celltype_averagecluster$RNA,"pbmc1_1_Celltype_averagecluster.txt",sep = "\t")
remove(pbmc1_1_Celltype_averagecluster,pbmc1_Celltype_averagecluster)

Idents(pbmc1)<-pbmc1@meta.data$Subset
DimPlot(pbmc1, reduction = "umap",label = TRUE)
DimPlot(pbmc1, reduction = "tsne",label = TRUE)

Idents(pbmc1_1)<-pbmc1_1@meta.data$Subset
DimPlot(pbmc1_1, reduction = "umap",label = TRUE)
DimPlot(pbmc1, reduction = "tsne",label = TRUE)

pbmc1_Subset.markers<-FindAllMarkers(pbmc1,only.pos = TRUE,logfc.threshold = 0,min.pct = 0)
pbmc1_1_Subset.markers<-FindAllMarkers(pbmc1_1,only.pos = TRUE,logfc.threshold = 0,min.pct = 0)
write.table(pbmc1_Subset.markers,"pbmc1_Subset.markers.txt",sep = "\t")
write.table(pbmc1_1_Subset.markers,"pbmc1_1_Subset.markers.txt",sep = "\t")
remove(pbmc1_1_Subset.markers,pbmc1_Subset.markers)

pbmc1_Subset_averagecluster<-AverageExpression(pbmc1,assays = "RNA")
pbmc1_1_Subset_averagecluster<-AverageExpression(pbmc1_1,assays = "RNA")
write.table(pbmc1_Subset_averagecluster$RNA,"pbmc1_Subset_averagecluster.txt",sep = "\t")
write.table(pbmc1_1_Subset_averagecluster$RNA,"pbmc1_1_Subset_averagecluster.txt",sep = "\t")
remove(pbmc1_1_Subset_averagecluster,pbmc1_Subset_averagecluster)

preGC_PB_GC1<-subset(pbmc1,idents = c("FCRL2/3high GC","preGC","LZ GC","prePB","GC","DZ GC","Cycling B","Plasmablast"))
GC1<-subset(pbmc1,idents = c("FCRL2/3high GC","LZ GC","prePB","GC","DZ GC","Cycling B"))

GC1.count<-GC1@assays$RNA@counts
GC1.metaeata<-GC1@meta.data
GC1_new<-CreateSeuratObject(counts = GC1.count, project = "pbmc3k", min.cells = 3, min.features = 200,meta.data = GC1.metaeata)
remove(GC1.count,GC1.metaeata)
GC1_new[["percent.mt"]] <- PercentageFeatureSet(GC1_new, pattern = "^MT-")
VlnPlot(GC1_new, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
GC1_new <- subset(GC1_new, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)
GC1_new <- NormalizeData(GC1_new)
GC1_new <- FindVariableFeatures(GC1_new, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(GC1_new)
GC1_new <- ScaleData(GC1_new, features = all.genes)
GC1_new <- RunPCA(GC1_new, features = VariableFeatures(object = GC1_new))
ElbowPlot(GC1_new,ndims = 50)
GC1_new <- FindNeighbors(GC1_new, dims = 1:20)
GC1_new <- FindClusters(GC1_new, resolution = 0.8)
GC1_new <- RunUMAP(GC1_new, dims = 1:20)
GC1_new <- RunTSNE(GC1_new, dims = 1:20)


Idents(GC1_new)<-GC1_new@meta.data$Subset
VlnPlot(GC1_new,features = c("MYC","IRF4","KDM6B","LY75"), pt.size = 0, ncol = 2)
DimPlot(GC1_new)
DimPlot(GC1_new,reduction = "tsne")
FeaturePlot(GC1_new,features = c("MYC","IRF4","KDM6B","LY75"),label = TRUE)
FeaturePlot(GC1_new,features = c("MYC","IRF4","KDM6B","LY75"),label = TRUE,reduction = "tsne")
GC1_new_Subset.markers<-FindAllMarkers(GC1_new,only.pos = TRUE,logfc.threshold = 0,min.pct = 0)
GC1_new_Subset_averagecluster<-AverageExpression(GC1_new,assays = "RNA")
write.table(GC1_new_Subset_averagecluster$RNA,"GC1_new_Subset_averagecluster.txt",sep = "\t")
write.table(GC1_new_Subset.markers,"GC1_new_Subset.markers.txt",sep = "\t")
remove(all.genes,GC1_new_Subset.markers,GC1_new_Subset_averagecluster)
remove(GC1,GC1_new,pbmc1,pbmc1_1,preGC_PB_GC1)

Idents(pbmc)<-pbmc@meta.data$seurat_clusters
Idents(pbmc_1)<-pbmc@meta.data$seurat_clusters

pbmc.markers<-FindAllMarkers(pbmc,only.pos = TRUE,logfc.threshold = 0,min.pct = 0)
pbmc_1.markers<-FindAllMarkers(pbmc_1,only.pos = TRUE,logfc.threshold = 0,min.pct = 0)
