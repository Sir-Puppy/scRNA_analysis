
# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(AnnotationHub)
library(ensembldb)

#read data for sample
GSM2883182_PLN_ <- read.csv("Data/GSM2883182_PLN++.csv")
GSM2883183_PLNr9c <- read.csv("Data/GSM2883183_PLNr9c.csv")

counts1 <- GSM2883182_PLN_
counts2 <- GSM2883183_PLNr9c

# Turn count matrix into a Seurat object (output is a Seurat object)
PLN_ <- CreateSeuratObject(GSM2883182_PLN_, min.features = 100)

PLNr9c <- CreateSeuratObject(GSM2883183_PLNr9c, min.features = 100)



# Create a merged Seurat object
merged_seurat <- merge(x = PLN_, 
                       y = PLNr9c, 
                       add.cell.id = c("PLN_", "PLNr9c"))

# Explore merged metadata
View(merged_seurat@meta.data)

# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

######################################################################################################################################
#We need to manually find the mito ratio. Human is likely there

# Connect to AnnotationHub
ah <- AnnotationHub()

# Access the Ensembl database for Homo sapiens and Mus musculus organisms
ahHS <- query(ah, 
              pattern = c("Homo sapiens", "EnsDb"), 
              ignore.case = TRUE)
ahMM <- query(ah, 
              pattern = c("Mus musculus", "EnsDb"), 
              ignore.case = TRUE)

# Check versions of databases available
ahHS %>% 
  mcols()
ahMM %>% 
  mcols()

# Acquire the latest annotation files for Homo sapiens and Mus musculus
idHS <- ahHS %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)
idMM <- ahMM %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

# Download the appropriate Ensembldb databases for Homo sapiens and Mus musculus
edbHS <- ah[[idHS]]
edbMM <- ah[[idMM]]

# Extract gene-level information from databases for Homo sapiens and Mus musculus
annotationsHS <- genes(edbHS, 
                       return.type = "data.frame")
annotationsMM <- genes(edbMM, 
                       return.type = "data.frame")    

# Select annotations of interest for Homo sapiens and Mus musculus
annotationsHS <- annotationsHS %>%
  dplyr::select(gene_id, gene_name, gene_biotype, seq_name, description, entrezid)
annotationsMM <- annotationsMM %>%
  dplyr::select(gene_id, gene_name, gene_biotype, seq_name, description, entrezid)

# Extract IDs for mitochondrial genes for Homo sapiens and Mus musculus
mtHS <- annotationsHS %>%
  dplyr::filter(seq_name == "MT") %>%
  dplyr::pull(gene_id)
mtMM <- annotationsMM %>%
  dplyr::filter(seq_name == "MT") %>%
  dplyr::pull(gene_id)

# Number of UMIs assigned to mitochondrial genes for Homo sapiens and Mus musculus
metadata$mtUMI_HS <- Matrix::colSums(counts1[which(rownames(counts1) %in% mtHS),], na.rm = T)
metadata$mtUMI_MM <- Matrix::colSums(counts2[which(rownames(counts2) %in% mtMM),], na.rm = T)

# Calculate of mitoRatio per cell for Homo sapiens and Mus musculus
metadata$mitoRatio_HS <- metadata$mtUMI_HS/metadata$nUMI
metadata$mitoRatio_MM <- metadata$mtUMI_MM/metadata$nUMI



########################################################################################################################################
#Strangely, no mitochondrial sequences were found, so we resorted to the following 3 lines of code

# Compute percent mito ratio
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

# Create metadata dataframe
metadata <- merged_seurat@meta.data

#######################################################################################################
#continued as usual from here

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^PLN_"))] <- "PLN_"
metadata$sample[which(str_detect(metadata$cells, "^PLNr9c"))] <- "PLNr9c"

# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata

# Create .RData object to load at any time
save(merged_seurat, file="data/merged_filtered_seurat.RData")

# Visualize the number of cell counts per sample
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)


# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)
