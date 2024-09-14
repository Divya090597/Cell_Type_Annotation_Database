

# Load Required Libraries -------------------------------------------------

library(SingleR)
library(scGate)
library(dplyr)
library(remotes)
library(UCell)
library(Seurat)
library(ggplot2)
library(purrr)

# Input Seurat object file ----------------------------------------------

seu_obj = readRDS("~/data/Kidney/ccRCC_Updated_Seurat_object_Cellcycle_score_SI_18854_Tumor.rds")


# USER INPUT INFORMATION --------------------------------------------------

# Input directory and project name

dir = "~/data/Kidney/GSE159115/Scgate result" 
Prj = "GSE159115"

# Output directories

output_dir = getwd()
meta_dir = "GSE159115_meta_data_combined"
cluster_annotation_dir = "GSE159115_cluster_annotation_percentage"


# scGate Annotation -------------------------------------------------------

# List of model names and corresponding scGate models

scGate_models_DB <-  get_scGateDB()
models <- list(
  CD4_TIL = scGate_models_DB$human$CD4_TIL,
  CD8_TIL = scGate_models_DB$human$CD8_TIL,
  generic = scGate_models_DB$human$generic,
  HiTME = scGate_models_DB$human$HiTME,
  TME_broad = scGate_models_DB$human$TME_broad,
  TME_HiRes = scGate_models_DB$human$TME_HiRes
)

# Create a list to store new columns
new_columns_list <- list()

# Loop through each model

for (model_name in names(models)) {
  # Get the current model
  my_scGate_model <- models[[model_name]]
  
  # Copy the original Seurat object to avoid modification
  seu_obj_copy <- seu_obj
  
  # Apply scGate to the copied Seurat object
  seu_obj_processed <- scGate(seu_obj_copy, model = my_scGate_model, pos.thr = 0.5, neg.thr = 0.5, ncores = 8)
  
  # Assign the processed Seurat object to a variable named according to the model name
  assign(paste0("seu_obj_", model_name), seu_obj_processed)
  
  # Extract metadata
  meta_data <- seu_obj_processed@meta.data
  
  # Identify columns that start with "is.pure_"
  pure_columns <- grep("^is.pure_", colnames(meta_data), value = TRUE)
  
  # Initialize the new column to store the features
  new_col_name <- paste0("scGate_", model_name, "_Celltype")
  meta_data[[new_col_name]] <- NA
  
  # Loop through each "is.pure_" column and update the new column
  for (col in pure_columns) {
    feature <- sub("^is.pure_", "", col)
    pure_barcodes <- rownames(meta_data)[meta_data[[col]] == "Pure"]
    for (barcode in pure_barcodes) {
      if (is.na(meta_data[barcode, new_col_name])) {
        meta_data[barcode, new_col_name] <- feature
      } else {
        meta_data[barcode, new_col_name] <- paste(meta_data[barcode, new_col_name], feature, sep = " ")
      }
    }
  }
  
  # Store the new column in the list
  new_columns_list[[new_col_name]] <- meta_data[, new_col_name, drop = F]
  
  # Save the metadata to a CSV file
  write.csv(meta_data, file = paste0(dir, "/meta_data_", model_name, ".csv"))
}


# Combine the original metadata with the new columns
meta_data_combined <- seu_obj@meta.data

for (new_col in new_columns_list) {
  meta_data_combined <- merge(meta_data_combined, new_col, by = "row.names", all.x = TRUE)
  rownames(meta_data_combined) <- meta_data_combined$Row.names
  meta_data_combined$Row.names <- NULL
}

# Update the original Seurat object with the combined metadata
seu_obj@meta.data <- meta_data_combined


# Automated Cluster Annotation using SingleR ------------------------------

# Use SingleR for automated annotation

reference_HPCA <- celldex::HumanPrimaryCellAtlasData()
singleR_result_HPCA <- SingleR(test = seu_obj@assays$RNA$data, ref = reference_HPCA, labels = reference_HPCA$label.main)
seu_obj$SingleR.labels.HPCA <- singleR_result_HPCA$labels

# Plot SingleR annotations

Idents(seu_obj) = "SingleR.labels.HPCA"
DimPlot(seu_obj, reduction = "umap") + ggtitle("SingleR annotation HPCA")

DimPlot(seu_obj_processed,group.by = c("CellOntology_name","Immune_UCell", "Lymphoid_UCell", "PanBcell_UCell" ,"Bcell_UCell", "APC_UCell") ,label = T, label.size = 3) 


# Calculate Percentage of Each Cell Type in Clusters  ---------------------


# Define columns with cell type information

cell_type_columns <- c(
  "scGate_generic_Celltype", "scGate_HiTME_Celltype",  
  "scGate_TME_broad_Celltype", "scGate_TME_HiRes_Celltype", "sctype_classification", "SingleR.labels.HPCA",  
  "scGate_CD4_TIL_Celltype", "scGate_CD8_TIL_Celltype"
)

# Store cell type percentages per cluster

cluster_percentage_list <- list()

# Loop through each cell type column

for (cell_type_col in cell_type_columns) {
  # Generate a df with cluster and cell type information
  cluster_celltype_df <- data.frame(
    #cluster = Idents(seu_obj), 
    cluster = seu_obj[['RNA_snn_res.0.5']],
    cell_type = seu_obj[[cell_type_col]]
  )
  cluster_celltype_df
  
  # Rename the cell_type column if it's not named correctly
  colnames(cluster_celltype_df)[2] <- "cell_type"
  
  # Convert NA values to a string so they are treated as a category
  cluster_celltype_df$cell_type <- ifelse(is.na(cluster_celltype_df$cell_type), "NA", cluster_celltype_df$cell_type)
  
  # Generate a table of cluster vs cell type counts 
  cluster_celltype_count <- table(cluster_celltype_df$RNA_snn_res.0.5, cluster_celltype_df$cell_type)
  
  # Convert counts to percentages 
  cluster_celltype_percentage <- prop.table(cluster_celltype_count, margin = 1) * 100
  
  # Convert to df for easier viewing
  cluster_celltype_df <- as.data.frame(cluster_celltype_percentage)
  colnames(cluster_celltype_df) <- c('cluster', 'cell_type', 'percentage')
  
  # Cell type percentage grouped by cluster
  cluster_percentage_strings <- cluster_celltype_df %>% 
    group_by(cluster) %>% 
    arrange(desc(percentage)) %>% 
    summarise(!!paste0(cell_type_col, "percentage") := paste0(cell_type, "(", round(percentage, 2), "%)", collapse = " "))
  
  # Store the result in the list
  cluster_percentage_list[[cell_type_col]] <- cluster_percentage_strings
}

# Combine all cluster percentage data into a single dataframe
cluster_percentage_strings_combined <- reduce(cluster_percentage_list, full_join, by = "cluster")

# View the final combined data frame
View(cluster_percentage_strings_combined)


# Save cluster-level annotation and metadata ------------------------------

# Generate cluster level annotation
write.csv(cluster_percentage_strings_combined, file = paste0(dir, "/", cluster_annotation_dir, ".csv"))

# Generate the combined metadata to a csv file
write.csv(seu_obj@meta.data, file = paste0(dir, "/", meta_dir, ".csv"))

saveRDS(seu_obj,"~/data/Kidney/Updated_ccRCC_Seu_Obj_SI_18854_Sctype_Scgate_SingleR.rds")


