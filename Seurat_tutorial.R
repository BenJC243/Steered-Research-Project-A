library(dplyr)
library(Seurat)
library(patchwork)

#loading the dataset]
pbmc.data <- Read10X(data.dir = "/home/solomon/filtered_gene_bc_matrices/hg19/")

#Initialise the Seurat Object with raw (non-normalised data)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

#Since most values in an scRNA-seq matrix are 0, Seurat uses a sparse-matrix 
#representation whenever possible. This results in significant memory and speed
#savings for Drop-seq/inDrop/10x data.
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

# The [[ operator can add columns to object metadata. 
#This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, 
#but can be used for anything calculated by the object, i.e. columns in object 
#metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#This is the QC line of code 
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#Can do this instead 
pbmc <- NormalizeData(pbmc) 

#Identifying highly variable features - feature selection
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot3 <- VariableFeaturePlot(pbmc)
plot4 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot3 + plot4


#Scaling the data 
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
#this method is longer than just doing pbmc <- ScaleData(pbmc) but it will prevent
#us from working with heatmaps in the future 

#Perform linear dimensional reduction, performing a PCA on scaled data 
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

#Visualising PCA results
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

#To overcome the extensive technical noise in  scRNA-seq data, Seurat clusters cells 
#based on their PCA scores, with each PC essentially representing a ‘metafeature’ 
#The top principal components therefore represent a robust compression of the dataset.
#However, how many components should we choose to include? 10? 20? 100?
 
#resampling test inspired by the JackStraw procedure. We randomly permute a subset of the data (1% by default)
#and rerun PCA, constructing a ‘null distribution’ of feature scores, and repeat this procedure. 
#We identify ‘significant’ PCs as those who have a strong enrichment of low p-value features.

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

JackStrawPlot(pbmc, dims = 1:15)
#Significant’ PCs will show a strong enrichment of features with low p-values (solid curve above the dashed line

#Cluster the cells
#construct a KNN graph based on the euclidean distance in PCA space
pbmc <- FindNeighbors(pbmc, dims = 1:10)
#apply modularity optimization techniques such as the Louvain algorithm (default) 
#iteratively group cells together, with the goal of optimizing the standard modularity function. 
pbmc <- FindClusters(pbmc, resolution = 0.5)

#Clusters can be found using the Idents() function
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)


#Run- non-linear dimensional reduction (UMAP/tSNE)
# If you haven't installed UMAP, you can do so via 
#reticulate::py_install(packages = 'umap-learn')

pbmc <- RunUMAP(pbmc, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap", label = TRUE)

saveRDS(pbmc, file = "/home/solomon/rstudio_objects/pbmc_tutorial.rds")

#Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
#Visualising marker expression, shows expression probability distrivbution accross the clusters
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))


pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

#we can use canonical markers to easily match the unbiased clustering to the known cell types
#check website
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(pbmc, file = "/home/solomon/rstudio_objects/pbmc_completed_tutorial.rds")
