library(R.utils)
library(Seurat)
library(dplyr)
library(Matrix)
library(patchwork)
library(tidyverse)
library(HGNChelper)
?BiocManager::install("rhdf5")
library(data.table)

#Create data folder
dir.create("data")

# URL of the file to be downloaded
file_url <- "ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174401/suppl/GSE174401%5Ffiltered%5Ffeature%5Fbc%5Fmatrix.h5"

# Destination file path
dest_file <- "data/GSE174401_filtered_feature_bc_matrix.h5"

# Download the file
download.file(url = file_url, destfile = dest_file, mode = "wb")

# Print a message to indicate the download is complete
cat("File downloaded to", dest_file, "\n")

# Decompress the file
gunzip("data/GSE174401_filtered_feature_bc_matrix.h5", remove = TRUE)


h5_file <- "data/GSE174401_filtered_feature_bc_matrix.h5"
Braindata <- Read10X_h5(h5_file)
Brain_data <- CreateSeuratObject(counts = Braindata, project = "Brain", min.cells = 3, min.features = 200)

#NORMALIZING THE DATA

Brain_data <- NormalizeData(Brain_data, normalization.method = "LogNormalize", scale.factor = 10000)

Brain_data <- FindVariableFeatures(Brain_data, selection.method = "vst", nfeatures = 2000)

# Data Scaling
Brain_data <- ScaleData(Brain_data)

# Perform PCA
Brain_data <- RunPCA(Brain_data)
Brain_data <- FindNeighbors(Brain_data, dims = 1:10)
Brain_data <- FindClusters(Brain_data, resolution = 0.3)
Brain_data <- RunUMAP(Brain_data, dims = 1:10)
DimPlot(Brain_data, reduction = "umap", label = TRUE, repel = TRUE)

saveRDS(Brain_data,"~/data/Brain_Met_Seurat.rds" )

#Find allmarkers
Markers = FindAllMarkers(Brain_data, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)

fwrite(Markers,"~/data/Brain_Met_Markers.csv")
library(openxlsx)
install.packages("openxlsx")

Sctype_DB <- read.xlsx("~/ScTypeDB_full.xlsx")
rm(gene_sets_prepare)
gs_list <- gene_sets_prepare(Sctype_DB, tissue)

