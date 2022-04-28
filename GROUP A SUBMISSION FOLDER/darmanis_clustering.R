install.packages("Rtsne")
install.packages("factoextra")
install.packages("mclust")


####Method 1#####
################
library(Rtsne)
library(mclust)
library(factoextra)

#Importing read count matrix 
braindata <- read.csv("darmanis_matrix.csv") #read count matrix 
rownames(braindata) <- braindata[,1] # assign rownames as column 1
braindata <- braindata[,-1] # remove column 1

set.seed(42)

# # filter out low-gene cells (often empty wells)
# braindata <- braindata[, colSums(braindata>0)>1.8e3]
# # remove genes that don't have many reads
# braindata <- braindata[rowSums(braindata)>10, ]
# # remove genes that are not seen in a sufficient number of cells
# braindata <- braindata[rowSums(braindata>0)>5, ]


#removing duplicate values from the data frame and turn it into a matrix object 
matrix_unique <-  as.matrix(unique(braindata))

#Creating a normalisation matrix for Rtsne
matrix_normalised <-  normalize_input(matrix_unique)

#transposing matrix, prevents genes being plotted against each other 
transpose_matrix <- t(matrix_normalised)
tsne_out <- Rtsne(as.matrix(dist(normalize_input(transpose_matrix))), theta=0.0)
tsne_y <- as.data.frame(tsne_out$Y)

set.seed(42)

# # filter out low-gene cells (often empty wells)
# braindata <- braindata[, colSums(braindata>0)>1.8e3]
# # remove genes that don't have many reads
# braindata <- braindata[rowSums(braindata)>10, ]
# # remove genes that are not seen in a sufficient number of cells
# braindata <- braindata[rowSums(braindata>0)>5, ]

#Creates plots
fit <- Mclust(tsne_y)
summary(fit)
plot(fit,what="BIC") #OPTIMAL CLUSTERS
plot(fit, what= "classification")#the higher the curve the better - VEV, CLUSTER

plot(fit, what ="BIC", xlab = "Number of Components")


# #Extra code 
# matrix_unique <-  as.matrix(unique(braindata))
# 
# #Creating a normalisation matrix for Rtsne
# matrix_normalised <-  normalize_input(matrix_unique)
# 
# #transposing matrix, prevents genes being plotted against each other 
# transpose_matrix <- t(matrix_normalised)
# 
# #Creating distance matrix of normalised data - returns 'dist' object 
# distance_matrix <- dist(as.matrix(transpose_matrix, method = "euclidean"))
# set.seed(42) #setting a seed allows for reproducible results
# tsne.model.1 <- Rtsne(distance_matrix, verbose=1,dims=3, num_threads=0, is_distance = TRUE, theta=0.0)
# d_tsne_1 = as.data.frame((tsne.model.1$Y)) #reduce dimensiosn held in a list
# 
# 
# 
# #Displaying a preferential model along with the number of clusters
# fit <- Mclust(d_tsne_1)#,G=9) #use G= to force the number of clusters
# summary(fit) #model VEV , 5 components
# plot(fit)
# 1
# plot(fit,what="BIC")
# plot(fit, what= "classification")#the higher the curve the better - VEV
# fviz_mclust_bic(fit)
# fviz_mclust(fit, 'classification')
# #fviz_cluster(fit, braindata[,-1],G=10)
