#Seurat clustering method

#install.packages("devtools")
#devtools::install_github("VCCRI/CIDR")
## install CIDR from https://github.com/VCCRI/CIDR
library(cidr)
library(mclust)

## Read in data 
## Tag tables were downloaded from the data repository NCBI Gene Expression Omnibus (GSE67835)
## Tag tables were combined into one table and hybrid cells have been excluded.
setwd("/home/solomon/Research_project")
brainTags <- read.csv("brainTags.csv")
rownames(brainTags) <- brainTags[,1]
brainTags <- brainTags[,-1]

## Read in annotation
info <- read.csv("SraRunTable.txt",sep="\t")
cellType <- info$cell_type_s[match(colnames(brainTags),info$Sample_Name_s)]
cellType <- factor(cellType)
types <- levels(cellType)

## Assign each cell type a color
scols <- c("red","blue","green","brown","pink","purple","darkgreen","grey")
cols <- rep(NA,length(cellType))
true_label <- rep(NA,length(cellType))
for (i in 1:length(cols)){
  cols[i] <- scols[which(types==cellType[i])]
  true_label[i] <- which(types==cellType[i])
}

## Clearer cell type names
types[3] <- "fetal quiescent neurons"
types[4] <- "fetal replicating neurons"
types[8] <- "oligodendrocyte precursor cells"

##  Standard principal component analysis using prcomp
priorTPM <- 1
brain10 <- brainTags[rowSums(brainTags)>10,]
brain10_lcpm <- log2(t(t(brain10)/colSums(brain10))*1000000+priorTPM)
pca <- prcomp(t(brain10_lcpm))
plot.new()
legend("center", legend=types, col=scols, pch=1)
plot(pca$x[,c(1,2)], col=cols, pch=1,
     xlab="PC1", ylab="PC2",
     main="Principal Component Analysis (prcomp)")


#### Seurat ####
################


BiocManager::install("Seurat") 

library(Seurat)
library(dplyr)
library(Matrix)
# Load the PBMC dataset

#pbmc.data <- read.csv("data_sc3_seurat/A_100.csv", row.names = 1);

start_time <- Sys.time()

# Examine the memory savings between regular and sparse matrices
dense.size <- object.size(x = as.matrix(x = brainTags))
sparse.size <- object.size(x = brainTags)

mincell=0 # Keep all genes expressed in >= mincell cells (parameter)
mingene=0 # Keep all cells with at least mingene detected genes (parameter)

pbmc <- CreateSeuratObject(brainTags)
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 100)
pbmc <- FindVariableFeatures(pbmc)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
pbmc <- RunPCA(pbmc, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
j=1; # a tunable parameter
pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:75, nn.eps = 0.5)
results <- FindClusters(object = pbmc, resolution = j)

end_time <- Sys.time()
t_Seurat = end_time - start_time
print(t_Seurat)

ARI_Seurat <- adjustedRandIndex(results@active.ident,cols)
ARI_Seurat
