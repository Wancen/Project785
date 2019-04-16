# Author: Weifang Liu
# Date: Apr 15

# Followed the Seurant pbmc3k tutorial: https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html

# Seurat on the e14 batch corrected data set 
library(dplyr)
library(Seurat)
library(ggplot2)
# e14 <- read.csv("/Users/Lugia/Documents/2019_Spring/Bios785/785project/sc_batch_corrected_data.csv")

# Load data 
load("/Users/Lugia/Documents/2019_Spring/Bios785/785project/e14.rda")
View(e14[1:5,1:5])
# Wwitch row names to genes 
row.names(e14) <- e14$X
# Delete the first column which are genes
e14 <- e14[,-1]
# Create a Seurat object 
e14 <- CreateSeuratObject(counts = e14, project = "mouse_cortex_e14")
# Check the Seurat object: 12982 features across 10931 samples
e14
dim(e14)

# Since the data is already filtered and median-normalized, we skip the standard pre-processing workflow and normalizing the data steps

# Feature selection: select the 2000 most variable genes for downstream analysis
e14_2000 <- FindVariableFeatures(object = e14, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = e14_2000), 10)
top10
# Plot variable features without labels
plot1 <- VariableFeaturePlot(object = e14_2000)
plot1

# Scale the data for PCA 
all.genes <- rownames(x=e14_2000)
e14_2000 <- ScaleData(object = e14_2000, features = all.genes)

# PCA
e14_2000 <- RunPCA(object = e14_2000, features = VariableFeatures(object = e14_2000))

# Examine and visualize PCA results a few different ways
print(x = e14_2000[["pca"]], dims = 1:5, nfeatures = 5)
# Visualize the first two PCs loadings
pdf("e14_2000_pca_loading.pdf")
VizDimLoadings(object = e14_2000, dims = 1:2, reduction = "pca")
dev.off()
# Plot PCA
pdf("PCA_plot.pdf")
DimPlot(object = e14_2000, reduction = "pca")
dev.off()
# Heatmaps
pdf("heatmap_p1")
DimHeatmap(object = e14_2000, dims = 1, cells = 500, balanced = TRUE)
dev.off()
pdf("heatmap_p1-9")
DimHeatmap(object = e14_2000, dims = 1:9, cells = 500, balanced = TRUE)
dev.off()
pdf("heatmap_p10-18")
DimHeatmap(object = e14_2000, dims = 10:21, cells = 500, balanced = TRUE)
dev.off()

# Determine the dimensionality of the dataset
e14_2000 <- JackStraw(object = e14_2000, num.replicate = 100, dim = 50)

e14_2000_score <- ScoreJackStraw(object = e14_2000, dims = 1:50)
pdf("JackStrawPlot")
JackStrawPlot(object = e14_2000_score, dims = 15:30)
dev.off()
pdf("ElbowPlot")
ElbowPlot(object = e14_2000, ndims = 50)
dev.off()

# try different number of PCs
# choose 30 for now according to the elbow plot

# Cluster the cells 
e14_cluster <- FindNeighbors(object = e14_2000_score, dims = 1:30)
e14_cluster <- FindClusters(object = e14_cluster, resolution = 0.5) # what is the resolution? 
# Look at cluster IDs of the first 5 cells
# Look at cluster IDs of the first 5 cells
head(x = Idents(object = e14_cluster), 5)

# Run non-linear dimensional reduction (UMAP/tSNE)
# Installed UMAP via reticulate::py_install(packages ='umap-learn')
e14_umap <- RunUMAP(object = e14_cluster, dims = 1:30)
# individual clusters
pdf("UMAP")
DimPlot(object = e14_umap, reduction = "umap")
dev.off()

# Finding differentially expressed features (cluster biomarkers)