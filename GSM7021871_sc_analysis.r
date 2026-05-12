# Single Cell Spatial Omics Analysis for GSM7021871
# Load necessary libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)

##### GSM7021871 Data Analysis #####
###### HCC Data Analysis ######

# Load the H5 file
hcc_counts <- Read10X_h5("./HCC/filtered_feature_bc_matrix.h5")
raw_hcc_counts <- Read10X_h5("./HCC/raw_feature_bc_matrix.h5")
# Create a standard Seurat Object
og_hcc_seurat <- CreateSeuratObject(counts = hcc_counts, project = "HCC_SingleCell")
raw_hcc_seurat <- CreateSeuratObject(counts = raw_hcc_counts, project = "HCC_SingleCell_Raw")

# Raw data analysis
## Quality Control
raw_hcc_seurat[["percent.mt"]] <- PercentageFeatureSet(raw_hcc_seurat, pattern = "^MT-")
raw_hcc_seurat[["percent.ribo"]] <- PercentageFeatureSet(raw_hcc_seurat, pattern = "^RPS|^RPL")

# Visualize QC metrics
qc_plot <- VlnPlot(raw_hcc_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)
print(qc_plot)