library(Seurat)
library(dplyr)
library(readr)
library(stringr)
library(Matrix)

# Load the CSV raw data
# Replace "path/to/your/data.csv" with the actual file path

#data1 <- read_csv("GSM2883182_PLN++.csv", col_types = cols()) 
data2 <- read_csv("GSM2883183_PLNr9c.csv", col_types = cols()) 

# Create a Seurat object
seurat_obj <- CreateSeuratObject(counts = data2)

# Perform quality control
# Filter cells based on minimum and maximum expressed genes and total counts
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 500 & nCount_RNA < 20000)

# Normalize data and identify highly variable genes
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Scale and center the data
seurat_obj <- ScaleData(seurat_obj)

# Perform dimensionality reduction using PCA
seurat_obj <- RunPCA(seurat_obj)

# Find neighbors and perform clustering using the Louvain algorithm

seurat_obj <- FindNeighbors(seurat_obj)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# TSNE visualization
seurat_obj <- RunTSNE(seurat_obj, tsne.method = "Rtsne", reduction = "pca")

# Plot the UMAP with cell clusters
DimPlot(seurat_obj, group.by = "seurat_clusters", reduction = "tsne")

# Access cluster information
cluster_ids <- seurat_obj$seurat_clusters
table(Idents(seurat_obj))

# Downstream analysis

# Save 
#saveRDS(seurat_obj, file = "GSM2883182_PLN++.rds")
