install.packages(c("SeuratObject","Seurat","cidr","dplyr","Matrix"))

library(SeuratObject)
library(Seurat)
library(cidr)
library(dplyr)
library(Matrix)

#Change this to your own directory 
#setwd("/home/solomon/Research_project")

darmanis <- read.csv("darmanis_matrix.csv")
rownames(darmanis) <- darmanis[,1]
braindata <- darmanis[,-1]



# Load the pbmc dataset

pbmc.data <- braindata

dim(braindata)



###Pre-processing workflow###

#Initalising the Seurat object
pbmc <- CreateSeuratObject(braindata)

dim(pbmc)

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Normalizing data
#pbmc <- NormalizeData(pbmc)
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 100)


#Identifying highly varaible features
pbmc <- FindVariableFeatures(pbmc)
#pbmc <- SCTransform(pbmc, assay= RNA, variable.features.n = 300)
top10 <- head(VariableFeatures(pbmc), 10) #10 most highly variable genes
#Plotting variable features with labels 
plot3 <- VariableFeaturePlot(pbmc)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
plot4


#Normalizing data
pbmc <- NormalizeData(pbmc)
#pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10

#Scaling the data
#pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
pbmc <- NormalizeData(pbmc)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
#pbmc <- ScaleData(pbmc)
pbmc <- NormalizeData(pbmc)
#Linear dimensional reduction -> PCA
#pbmc <- RunPCA(pbmc, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:10, nfeatures = 5)# Examine and visualize PCA results a few different ways

#Visualising cells and features that define the PCA
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1:10, cells = 500, balanced = TRUE)

dim(pbmc)

#Determinig dimensionality of the dataset basesd on PCA score
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

JackStrawPlot(pbmc, dims = 1:20)

#Clustering the Cells
pbmc <- FindNeighbors(pbmc, dims = 1:11)
pbmc <- FindClusters(pbmc, resolution = 0.5) 

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)


#Non-linear dimensional reduction using UMAP and tsne
pbmc <- RunUMAP(pbmc, dims = 1:11)
DimPlot(pbmc, reduction = "umap", label = TRUE)
# saveRDS(pbmc, 'saverds/pipeline1.rds')
# 
# pbmc <- RunTSNE(pbmc, dims=1:34)
# # tSNE 
# tSNE_data <- RunTSNE(cluster_data, dims=1:34)
# DimPlot(tSNE_data, reduction='tsne')

#Still need to fix this bit so we can select unbiased what cell type each cluster is
# find all markers of cluster 1
# cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, min.pct=-0.25)
# cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
# head(cluster0.markers, n = 5)
# cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct =0.25)
# cluster3.markers <- FindMarkers(pbmc, ident.1 = 3, min.pct =0.25)
# cluster4.markers <- FindMarkers(pbmc, ident.1 = 4, min.pct=-0.25)
# cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, min.pct=-0.25)
# cluster6.markers <- FindMarkers(pbmc, ident.1 = 6, min.pct=-0.25)
# 
# cluster0.markers
# 
# 
# # find all markers distinguishing cluster 5 from clusters 0 and 3
# Distinguishing.cluster0 <- FindMarkers(pbmc, ident.1 = 0, ident.2 = c(6, 3, 1, 2, 4, 5), min.pct = 0.25)
# head(Distinguishing.cluster0, n = 5)
# 
# Distinguishing.cluster1 <- FindMarkers(pbmc, ident.1 = 1, ident.2 = c(6, 3, 0, 2, 4, 5), min.pct = 0.25)
# head(Distinguishing.cluster1, n = 5)
# 
# Distinguishing.cluster2 <- FindMarkers(pbmc, ident.1 = 2, ident.2 = c(6, 3, 0, 1, 4, 5), min.pct = 0.25)
# head(Distinguishing.cluster2, n = 20)
# 
# Distinguishing.cluster3 <- FindMarkers(pbmc, ident.1 = 3, ident.2 = c(6, 2, 0, 1, 4, 5), min.pct = 0.25)
# head(Distinguishing.cluster3, n = 5)
# 
# Distinguishing.cluster4 <- FindMarkers(pbmc, ident.1 = 4, ident.2 = c(6, 2, 0, 1, 3, 5), min.pct = 0.25)
# head(Distinguishing.cluster4, n = 5)
# 
# Distinguishing.cluster5 <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(6, 2, 0, 1, 3, 4), min.pct = 0.25)
# head(Distinguishing.cluster5, n = 5)
# 
# Distinguishing.cluster6 <- FindMarkers(pbmc, ident.1 = 6, ident.2 = c(5, 2, 0, 1, 3, 4), min.pct = 0.25)
# head(Distinguishing.cluster6, n = 5)
# 
# saveRDS(pbmc ,'saverds/pipeline1.rds', compress = FALSE)
# 
# table(Idents(pbmc))

