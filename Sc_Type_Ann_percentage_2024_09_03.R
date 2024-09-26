
# Load required libraries -------------------------------------------------

library(openxlsx)
library(Seurat)
library(dplyr)
library(scales)

# load gene set preparation function --------------------------------------


gene_sets_prepare <- function(path_to_db_file, cell_type){
  
  cell_markers = openxlsx::read.xlsx(path_to_db_file)
  cell_markers = cell_markers[cell_markers$tissueType == cell_type,] 
  cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1); cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  # correct gene symbols from the given DB (up-genes)
  cell_markers$geneSymbolmore1 = sapply(1:nrow(cell_markers), function(i){
    
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore1[i],",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    
    if(length(markers_all) > 0){
      suppressMessages({markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))})
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })
  
  # correct gene symbols from the given DB (down-genes)
  cell_markers$geneSymbolmore2 = sapply(1:nrow(cell_markers), function(i){
    
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore2[i],",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    
    if(length(markers_all) > 0){
      suppressMessages({markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))})
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })
  
  cell_markers$geneSymbolmore1 = gsub("///",",",cell_markers$geneSymbolmore1);cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gsub("///",",",cell_markers$geneSymbolmore2);cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  gs = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]),",")))); names(gs) = cell_markers$cellName
  gs2 = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]),",")))); names(gs2) = cell_markers$cellName
  
  list(gs_positive = gs, gs_negative = gs2)
}


# load cell type annotation function --------------------------------------


sctype_score <- function(scRNAseqData, scaled = !0, gs, gs2 = NULL, gene_names_to_uppercase = !0, ...){
  
  # check input matrix
  if(!is.matrix(scRNAseqData)){
    warning("scRNAseqData doesn't seem to be a matrix")
  } else {
    if(sum(dim(scRNAseqData))==0){
      warning("The dimension of input scRNAseqData matrix equals to 0, is it an empty matrix?")
    }
  }
  
  # marker sensitivity
  marker_stat = sort(table(unlist(gs)), decreasing = T); 
  marker_sensitivity = data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_stat), to = c(0,1), from = c(length(gs),1)),
                                  gene_ = names(marker_stat), stringsAsFactors = !1)
  
  # convert gene names to Uppercase
  if(gene_names_to_uppercase){
    rownames(scRNAseqData) = toupper(rownames(scRNAseqData));
  }
  
  # subselect genes only found in data
  names_gs_cp = names(gs); names_gs_2_cp = names(gs2);
  gs = lapply(1:length(gs), function(d_){ 
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  gs2 = lapply(1:length(gs2), function(d_){ 
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs2[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  names(gs) = names_gs_cp; names(gs2) = names_gs_2_cp;
  cell_markers_genes_score = marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(gs)),]
  
  # z-scale if not
  if(!scaled) Z <- t(scale(t(scRNAseqData))) else Z <- scRNAseqData
  
  # multiple by marker sensitivity
  for(jj in 1:nrow(cell_markers_genes_score)){
    Z[cell_markers_genes_score[jj,"gene_"], ] = Z[cell_markers_genes_score[jj,"gene_"], ] * cell_markers_genes_score[jj, "score_marker_sensitivity"]
  }
  
  # subselect only with marker genes
  Z = Z[unique(c(unlist(gs),unlist(gs2))), ]
  
  # combine scores
  es = do.call("rbind", lapply(names(gs), function(gss_){ 
    sapply(1:ncol(Z), function(j) {
      gs_z = Z[gs[[gss_]], j]; gz_2 = Z[gs2[[gss_]], j] * -1
      sum_t1 = (sum(gs_z) / sqrt(length(gs_z))); sum_t2 = sum(gz_2) / sqrt(length(gz_2));
      if(is.na(sum_t2)){
        sum_t2 = 0;
      }
      sum_t1 + sum_t2
    })
  })) 
  
  dimnames(es) = list(names(gs), colnames(Z))
  es.max <- es[!apply(is.na(es) | es == "", 1, all),] # remove na rows
  
  es.max
}


#  User Input and Output --------------------------------------------------

# Input DB and study tissue (for sc_Type package)

db_ <- "ScTypeDB_full.xlsx";
tissue <- "Kidney" 

# Input Seurat object file ----------------------------------------------

seu_obj = readRDS("~/data/Kidney/ccRCC_Updated_Seurat_object_Cellcycle_score_SI_18854_Tumor.rds")


# USER INPUT INFORMATION --------------------------------------------------

# Input directory and project name

dir = "~/data/Kidney/GSE159115/Sc_Type_result" 
Prj = "GSE159115"

# Output directories

output_dir = getwd()
meta_dir = "GSE159115_meta_data_combined"
cluster_annotation_dir = "GSE159115_cluster_annotation_percentage"


# Input Seurat object file ----------------------------------------------

seu_obj = readRDS("~/data/Kidney/ccRCC_Updated_Seurat_object_Cellcycle_score_SI_18854_Tumor.rds")


# Prepare gene sets

gs_list <- gene_sets_prepare(db_, tissue)

db_ = openxlsx::read.xlsx("~/ScTypeDB_full.xlsx")

# Extract scaled data from Seurat object

scRNAseqData_scaled <- as.matrix(seu_obj[["RNA"]]$scale.data)

# Run ScType to score gene sets

es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# Merge scores by cluster

cL_resutls <- do.call("rbind", lapply(unique(seu_obj@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(seu_obj@meta.data[seu_obj@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seu_obj@meta.data$seurat_clusters==cl)), 10)
}))

# Group results by cluster and select top cell types

sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# Set low-confident (low ScType score) clusters to "unknown"

sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
print(sctype_scores[,1:3])

# Update Seurat object metadata

seu_obj@meta.data$sctype_classification = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  seu_obj@meta.data$sctype_classification[seu_obj@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

# Plot the results

DimPlot(seu_obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype_classification')        


#Calculate and Save Cell Type Percentages in Clusters -----------------


cell_type_columns <-  "sctype_classification"

Sc_type_cluster_percentage_list <- list()

generate_cluster_celltype_df <- function(seu_obj, cell_type_column) {
  # Generate a data frame with cluster and cell type information
  cluster_celltype_df <- data.frame(
    cluster = seu_obj$RNA_snn_res.0.5, 
    cell_type = seu_obj[[cell_type_column]]
  )
  
  # Return the generated data frame
  return(cluster_celltype_df)
}

cell_type_column <- "sctype_classification"
cluster_celltype_df <- generate_cluster_celltype_df(seu_obj, cell_type_column)

# Rename the cell_type column if it's not named correctly

colnames(cluster_celltype_df)[2] <- "cell_type"

# Convert NA values to a string so they are treated as a category

cluster_celltype_df$cell_type <- ifelse(is.na(cluster_celltype_df$cell_type), "NA", cluster_celltype_df$cell_type)

# Generate a table of cluster vs cell type counts 
cluster_celltype_count <- table(cluster_celltype_df$cluster, cluster_celltype_df$cell_type)

# Convert counts to percentages 
cluster_celltype_percentage <- prop.table(cluster_celltype_count, margin = 1) * 100

# Convert to df for easier viewing

cluster_celltype_df <- as.data.frame(cluster_celltype_percentage)
colnames(cluster_celltype_df) <- c('cluster', 'cell_type', 'percentage')

# Cell type percentage grouped by cluster

cluster_percentage_strings <- cluster_celltype_df %>% 
  group_by(cluster) %>% 
  arrange(desc(percentage)) %>% 
  summarise(!!paste0(cell_type_column, "_percentage") := paste0(cell_type, "(", round(percentage, 2), "%)", collapse = " "))

# Store the result in the list
Sc_type_cluster_percentage_list[[cell_type_column]] <- cluster_percentage_strings

# Generate cluster level annotation
write.csv(cluster_percentage_strings, file = paste0(dir, "/", cluster_annotation_dir, ".csv"))
