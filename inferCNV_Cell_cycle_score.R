
#  Load Required Libraries ------------------------------------------------

library(infercnv)
library(data.table)
library(dplyr)
library(biomaRt)
library(Seurat)
library(ggplot2)
library(tidyr)


# Load Seurat Objects -----------------------------------------------------

# Load Seurat objects for tumor and benign adjacent samples

seurat_obj_Tumor = readRDS("~/Kidney/ccRCC_seurat_object_Tumor.rds")

seurat_obj_BA = readRDS("~/Kidney/ccRCC_seurat_object_Benign_adjacent.rds")


#  Gene List Preparation Using BiomaRt ------------------------------------

# Get all gene names from Seurat object
all_genes_T <- rownames(seurat_obj_Tumor)

# Use Ensembl to retrieve gene information
human = useEnsembl("genes", dataset = "hsapiens_gene_ensembl")

# Generate a dataframe with gene information

gene_lists = getBM(attributes = c('hgnc_symbol','chromosome_name','start_position','end_position'), mart = human)

# Filter to keep genes only on main chromosomes

gene_lists = gene_lists %>%
  filter(chromosome_name %in% c(1:22, "X", "Y", "MT"),
         !(is.na(hgnc_symbol) | hgnc_symbol == ""))


#  Filter and Prepare Gene Order Data for inferCNV ------------------------


# Identify intersections and filter gene locations

gene_locs = gene_lists %>%
  filter(hgnc_symbol %in% all_genes_T) %>%
  group_by(hgnc_symbol) %>%
  mutate(entries = n()) %>%
  #because there are 2 genes that are on x and y, remove y
  filter(!c(entries == 2 & chromosome_name == "Y"))

#missing genes
all_genes_T[!c(all_genes_T %in% gene_locs_T$hgnc_symbol)]

#creating gene_order_file for infercnv

gene_locs <- gene_lists %>%
  filter(hgnc_symbol %in% all_genes_T) %>%
  group_by(hgnc_symbol) %>%
  mutate(entries = n()) %>%
  filter(!c(entries == 2 & chromosome_name == "Y")) %>%  # Remove entries on Y chromosome if there are duplicates
  ungroup() %>%
  dplyr::select(-entries) %>%
  distinct(hgnc_symbol, .keep_all = TRUE) %>%  # Keep only the first occurrence
  column_to_rownames("hgnc_symbol") %>%
  mutate(chromosome_name = paste0("chr", chromosome_name),
         chromosome_name = factor(chromosome_name, levels = paste0("chr", c(1:22, "X", "Y"))))

# Check for unique row names
anyDuplicated(rownames(gene_locs))


# Prepare Expression Matrix and Annotation Data  --------------------------


#creating annotation file, save without column names
gene_annotation = seurat_obj_Tumor@meta.data %>%
  dplyr::select(anno) 


#expression matrix of raw counts
express = seurat_obj_Tumor@assays$RNA$counts %>% 
  data.frame(check.names = FALSE)

# Create inferCNV Object ----------------------------------------------------------------

#create object
infercnv_object <- infercnv::CreateInfercnvObject(raw_counts_matrix=express, 
                                                  gene_order_file=gene_locs,
                                                  annotations_file=gene_annotation,
                                                  ref_group_names=NULL,
                                                  min_max_counts_per_cell = c(18, Inf))

express = infercnv_object@count.data
gene_locs = infercnv_object@gene_order
gene_annotation = gene_annotation %>%
  filter(row.names(.) %in% colnames(express))

#new object
infercnv_object <- infercnv::CreateInfercnvObject(raw_counts_matrix=express, 
                                                  gene_order_file=gene_locs,
                                                  annotations_file=gene_annotation,
                                                  ref_group_names=NULL,
                                                  min_max_counts_per_cell = c(18, Inf))

# perform infercnv operations to reveal cnv signal
infercnv_object = infercnv::run(infercnv_object,
                                cutoff=0.1, 
                                out_dir="~/data/Kidney/infercnv/Tumor",  # dir is auto-created for storing outputs
                                cluster_by_groups=T,   # cluster
                                denoise=T,
                                HMM=F #based on infercnv issue 234
)


#goes circular here

infercnv_object@reference_grouped_cell_indices #normal
infercnv_object@observation_grouped_cell_indices #tumor
Collapse

gene_order_file = write.csv(gene_locs,"~/data/Kidney/gene_order_file.csv")


# View the overall structure of the infercnv object
str(infercnv_object)

# List all the slots in the infercnv object
slotNames(infercnv_object)

# Access specific data within the object
head(infercnv_object@count.data)   # Access raw counts matrix
head(infercnv_object@gene_order)          # Access gene order information


# Fill missing chromosome names with "chrM" for mitochondrial genes
gene_locs$chr[is.na(gene_locs$chr)] <- "chrM"

# Convert the 'chr' column to character
gene_locs$chr <- as.character(gene_locs$chr)

# Replace missing values with "chrM"
gene_locs$chr[is.na(gene_locs$chr)] <- "chrM"

# (Optional) Convert back to factor if needed
gene_locs$chr <- factor(gene_locs$chr)

# Verify that there are no missing values remaining
colSums(is.na(gene_locs))

# Check the unique values in the 'chr' column
unique(gene_locs$chr)

# Rerun the analysis
infercnv_object = infercnv::run(infercnv_object,
                                cutoff=0.1, 
                                out_dir="~/data/Kidney/infercnv/Tumor",  # dir is auto-created for storing outputs
                                cluster_by_groups=T,   # cluster
                                denoise=T,
                                HMM=F #based on infercnv issue 234
)


# Calculate CNV Scores ----------------------------------------------------

cnv_matrix_Tumor <- infercnv_object@expr.data  # Extract the CNV matrix
cnv_score_Tumor <- apply(cnv_matrix_Tumor, 2, function(cell) sum(abs(cell - 2)))  # Calculate deviation score for each cell

View(cnv_score_Tumor)

cnv_score_Tumor = infercnv_object@expr.data %>% colMeans() %>% data.frame() %>% rename("CNV_Means" = 1) %>% rownames_to_column("cell_names")


metadata = seurat_obj_Tumor@meta.data
metadata$cell_names <- rownames(metadata)

# Merge CNV scores with metadata
merged_meta <- merge(metadata, cnv_score_Tumor, by = "cell_names")

seurat_obj_Tumor@meta.data = merged_meta_T

# Violin plot of CNV scores by cell type
seurat_obj_Tumor@meta.data %>% 
       rownames_to_column("cell_names") %>%
       full_join(infercnv_object@expr.data %>% colMeans() %>% data.frame() %>% rename("CNV_Means" = 1) %>% rownames_to_column("cell_names")) %>% 
       ggplot(aes(x = sctype_classification, y = CNV_Means, fill = sctype_classification)) +
       geom_hline(yintercept = 1)+
       geom_violin(trim = FALSE) +
       labs(title = "Violin Plot of CNV Scores", x = "Cell Type", y = "CNV Score") +
       theme_minimal()



# Calculate Cell Cycle score ----------------------------------------------

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seurat_obj_cellcycle <- CellCycleScoring(seurat_obj_Tumor, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seurat_obj_cellcycle@meta.data


# Plot UMAP colored by CNV mean scores

p1 <- FeaturePlot(seurat_obj, features = "CNV_Means") +
  scale_color_gradient(low = "yellow", high = "blue") +
  ggtitle("UMAP Plot of CNV Scores_Tumor")

# Plot UMAP colored by cell type classification
p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "sctype_classification") +
  ggtitle("UMAP Plot of Cell Type Classification")

p1 + p2
p1



FeaturePlot(seurat_obj_cellcycle, features = c("CNV_Means", "EPCAM", "MKI67"))

DimPlot(seurat_obj_cellcycle,reduction = "umap",  group.by = ("Phase"))

P1 = DimPlot(seurat_obj_cellcycle,reduction = "umap",  group.by = "sctype_classification")
  
P1

seurat_obj_cellcycle@meta.data

# Display the plot

DimPlot(seurat_obj,reduction = "umap",  group.by = "anno", label = T)

View(seurat_obj@meta.data)

# Save the updated Seurat objects and inferCNV object
saveRDS(seurat_obj,"~/data/Kidney/ccRCC_Updated_Seurat_Object_SI_18854_Tumor.rds")        
saveRDS(seurat_obj_cellcycle, "~/data/Kidney/ccRCC_Updated_Seurat_Object_Cellcycle_score_SI_18854_Tumor.rds")
saveRDS(infercnv_object_T,"~/data/infercnv/Tumor/ccRCC_InferCNV_Object_SI_18854_Tumor.rds")
