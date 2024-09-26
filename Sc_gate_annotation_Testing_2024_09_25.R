# Load Required Libraries

library(dplyr)
library(Seurat)
library(patchwork)
library(scGate)
library(plotly)
library(purrr)

# Data imported : "~/data/Kidney/ccRCC_seurat_object_Tumor.rds

seurat_obj = readRDS("~/data/Kidney/ccRCC_seurat_object_Tumor.rds")
my_scGate_model_gene_list1 <- scGate::load_scGate_model("Scgate_Immunecell_example_genelist_2024_09_25.tsv.txt")

obj <- scGate(seu_obj, model = my_scGate_model_gene_list1, save.levels = TRUE)
obj2 = scGate(seu_obj, model = scGate_models_DB$human$generic, save.levels = TRUE)
DimPlot(obj, cols = palette) + theme(aspect.ratio = 1)
plots <- scGate::plot_levels(seu_obj)
wrap_plots(plots, ncol = 2)

# My model

model = my_scGate_model_gene_list1
#models <- list(
  #Divya_model = my_scGate_model_gene_list1, 
  #Macrophage = scGate_models_DB$human$generic$Macrophage
#)

NewModel <- list(
  Divya_model = my_scGate_model_gene_list1
  #CD4Tcell= scGate_models_DB$human$generic$CD4T,
  #CD8Tcell= scGate_models_DB$human$generic$CD8T
  )

models <- list(
  NewModel = NewModel
)
# Create a list to store new columns
new_columns_list <- list()

# Loop through each model
model_name = "Divya_model"
for (model_name in names(models)) {
  # Get the current model
   my_scGate_model <- models[[model_name]]
  
  # Copy the original Seurat object to avoid modification
  seu_obj_copy <- seurat_obj
  
  # Apply scGate to the copied Seurat object
  seu_obj_processed <- scGate(seu_obj_copy, model = my_scGate_model, pos.thr = 0.2, neg.thr = 0.2, ncores = 4)
  
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
View(as.data.frame(new_columns_list))

# Combine the original metadata with the new columns
meta_data_combined <- seurat_obj@meta.data

for (new_col in new_columns_list) {
  meta_data_combined <- merge(meta_data_combined, new_col, by = "row.names", all.x = TRUE)
  rownames(meta_data_combined) <- meta_data_combined$Row.names
  meta_data_combined$Row.names <- NULL
}

# Update the original Seurat object with the combined metadata
seu_obj@meta.data <- meta_data_combined

