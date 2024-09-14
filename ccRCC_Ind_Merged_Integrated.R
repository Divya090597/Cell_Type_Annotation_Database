

# Load Required Libraries -------------------------------------------------


library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(rhdf5)


# GSE159115_SI_18854_ccRCC_Tumor_Tissue ------------------------------------------------


# Define the path to your .tar file
tar_file <- "~/GSE159115_RAW.tar"  

# Extract the .tar file to the current working directory
untar(tar_file)

# Define the path to the .h5 file
h5_file <- "~/GSM4819725_SI_18854_filtered_gene_bc_matrices_h5.h5"

# Load data using Read10X_h5
ccRCC_data <- Read10X_h5(h5_file)

# Create Seurat object for tumor tissue
Seurat_object_Tumor <- CreateSeuratObject(counts = ccRCC_data, project = "Renal", min.cells = 3, min.features = 200)

# Normalize the data
Seurat_object_Tumor <- NormalizeData(Seurat_object_Tumor, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify variable features
Seurat_object_Tumor <- FindVariableFeatures(Seurat_object_Tumor, selection.method = "vst", nfeatures = 2000)

# Scale the data
Seurat_object_Tumor <- ScaleData(Seurat_object_Tumor)

# Perform Principal Component Analysis (PCA)
Seurat_object_Tumor <- RunPCA(Seurat_object_Tumor)

# Find neighbors and clusters
Seurat_object_Tumor <- FindNeighbors(Seurat_object_Tumor, dims = 1:10)
Seurat_object_Tumor <- FindClusters(Seurat_object_Tumor, resolution = 0.5)

# Run UMAP for visualization
Seurat_object_Tumor <- RunUMAP(Seurat_object_Tumor, dims = 1:10)

# Plot UMAP with clusters
DimPlot(Seurat_object_Tumor, reduction = "umap", group.by = "RNA_snn_res.0.5", label = TRUE)



# GSE159115_SI_18855_ccRCC_Benign_Adjacent_Tissue ------------------------------------------------


# Define the path to the .h5 file
h5_file <- "~/GSM4819725_SI_18855_filtered_gene_bc_matrices_h5.h5"

# Load data using Read10X_h5
ccRCC_BA <- Read10X_h5(h5_file)

# Create Seurat object for tumor tissue
Seurat_object_BA <- CreateSeuratObject(counts = ccRCC_BA, project = "Renal", min.cells = 3, min.features = 200)

# Normalize the data
Seurat_object_TSeurat_object_BAumor <- NormalizeData(Seurat_object_BA, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify variable features
Seurat_object_BA <- FindVariableFeatures(Seurat_object_BA, selection.method = "vst", nfeatures = 2000)

# Scale the data
Seurat_object_BA <- ScaleData(Seurat_object_BA)

# Perform Principal Component Analysis (PCA)
Seurat_object_BA <- RunPCA(Seurat_object_BA)

# Find neighbors and clusters
Seurat_object_BA <- FindNeighbors(Seurat_object_BA, dims = 1:10)
Seurat_object_BA <- FindClusters(Seurat_object_BA, resolution = 0.5)

# Run UMAP for visualization
Seurat_object_BA <- RunUMAP(Seurat_object_BA, dims = 1:10)

# Plot UMAP with clusters
DimPlot(Seurat_object_BA, reduction = "umap", group.by = "RNA_snn_res.0.5", label = TRUE)



# Merge tumor and benign adjacent Seurat objects --------------------------


merged_seurat <- merge(Seurat_object_Tumor, y = Seurat_object_BA,add.cell.ids = c("Tumor", "Benign_Adjacent") ,project = "Merged")

# Save the merged Seurat object
saveRDS(merged_seurat, "~/Kidney/ccRCC_seurat_object_Tumor_Benign_adjacent_merged.rds")

# Normalize the merged data
merged_seurat <- NormalizeData(merged_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify variable features in the merged data
merged_seurat <- FindVariableFeatures(merged_seurat, selection.method = "vst", nfeatures = 2000)

# Scale the merged data
merged_seurat <- ScaleData(merged_seurat)

# Perform PCA on the merged data
merged_seurat <- RunPCA(merged_seurat)

# Find neighbors and clusters in the merged data
merged_seurat <- FindNeighbors(merged_seurat, dims = 1:10)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.5)

# Run UMAP for visualization of the merged data
merged_seurat <- RunUMAP(merged_seurat, dims = 1:10)

# Plot UMAP with clusters for the merged data
DimPlot(merged_seurat, reduction = "umap", group.by = "RNA_snn_res.0.5", label = TRUE)



# Define marker genes
genes <- c("IGKC", "IGHG1", "CD3E", "CD4", "CD8B", "CD79A", "MS4A1", "TPSB2", 
           "GZMB", "NKG7", "FCGR3A", "CD14", "CD19")

# Define genes to plot
genes_to_plot <- c("IGFBP3", "IGFBP5", "EPCAM", "VEGFA", "ANGPTL4", "PECAM1", 
                   "VCAM1", "ITGB8", "SLC17A3", "CD24", "BNIP3", "CRYAB", 
                   "MKI67", "PROM1", "CA9")

# Plot gene expression features
FeaturePlot(Urine_seurat, features = genes_to_plot, cols = c("lightgrey", "red"), pt.size = 1.5) +
  theme_minimal()

saveRDS(Kidney_seurat, "~/Data/Kidneycancer_SCRNASeq.rds")



# Integration -------------------------------------------------------------



#Split the object into a list of 2 repeats
mergeddata_list <- SplitObject(merged_seurat, split.by = 'orig.ident')

#normalize and identify variable features for each dataset independently
mergeddata_list <- lapply(X = mergeddata_list, FUN = function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = list(Tumor, Benign_Adjacent))

# Find integration anchors
integration_anchors <- FindIntegrationAnchors(object.list = list(Tumor, Benign_Adjacent), anchor.features = features)

# Perform integration to create an 'integrated' data assay
integrated_data <- IntegrateData(anchorset = integration_anchors, dims = 1:30)

# Perform a integration analysis
DefaultAssay(integrated_data) <- "integrated"

integrated_data <- ScaleData(integrated_data, verbose = FALSE)
integrated_data <- RunPCA(integrated_data, npcs = 30, verbose = FALSE)
integrated_data <- FindNeighbors(integrated_data, reduction = "pca", dims = 1:30)
integrated_data <- FindClusters(integrated_data, resolution = 0.5)
integrated_data <- RunUMAP(integrated_data, reduction = "pca", dims = 1:30)

# Visualization
DimPlot(integrated_data, reduction = "umap", label = TRUE)

# Compare
plot_1 <- DimPlot(mergeddata, reduction = "umap", group.by = 'orig.ident')
plot_2 <- DimPlot(integrated_data, reduction = "umap", group.by = 'orig.ident')
plot_1+plot_2







