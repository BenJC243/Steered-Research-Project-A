library(SeuratObject)
library(Seurat)
library(cidr)
library(dplyr)
library(Matrix)

#Change this to your own directory 
setwd("/home/solomon/Research_project")

#Reading in data
braindata <- read.csv("darmanis_matrix.csv") #read count matrix 
rownames(braindata) <- braindata[,1]
braindata <- braindata[,-1]

###Seurat Clustering####

# Load the PBMC dataset

pbmc.data <- braindata



#Initalising the Seurat object
pbmc <- CreateSeuratObject(braindata)

###Pre-processing workflow###

#Initalising the Seurat object
pbmc <- CreateSeuratObject(braindata)

#Quality control
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-") #mitochondrial QC metrics
#pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Normalizing data
pmbc <- NormalizeData(pbmc)
#pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 100)


#Identifying highly varaible features
pbmc <- FindVariableFeatures(pbmc)
top10 <- head(VariableFeatures(pbmc), 10) #10 most highly variable genes
#Plotting variable features with labels 
plot3 <- VariableFeaturePlot(pbmc)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
plot3 + plot4

#Normalizing data
pmbc <- NormalizeData(pbmc)
#pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10

#Scaling the data
#pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
pmbc <- NormalizeData(pbmc)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
#pbmc <- ScaleData(pbmc)

#Linear dimensional reduction -> PCA
#pbmc <- RunPCA(pbmc, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
pbmc <- RunPCA(pmbc, features = VariableFeatures(object = pmbc))
print(pbmc[["pca"]], dims = 1:10, nfeatures = 5)# Examine and visualize PCA results a few different ways

#Visualising cells and features that define the PCA
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1:10, cells = 500, balanced = TRUE)


#Determinig dimensionality of the dataset basesd on PCA score
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

JackStrawPlot(pbmc, dims = 1:15)

#Clustering the Cells
pbmc <- FindNeighbors(pbmc, dims = 1:12)
pbmc <- FindClusters(pbmc, resolution = 0.5) 

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)


#Non-linear dimensional reduction using UMAP and tsne
pbmc <- RunUMAP(pbmc, dims = 1:12)
DimPlot(pbmc, reduction = "umap", label = TRUE)

#Still need to fix this bit so we can select unbiased what cell type each cluster is
# find all markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)




