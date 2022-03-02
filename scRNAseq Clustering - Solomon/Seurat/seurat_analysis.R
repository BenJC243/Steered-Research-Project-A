#Seurat clustering method

#install.packages("devtools")
#devtools::install_github("VCCRI/CIDR")
## install CIDR from https://github.com/VCCRI/CIDR
library(cidr)

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

