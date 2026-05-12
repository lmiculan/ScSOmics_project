# Single Cell Spatial Omics Analysis for GSM7021871
# Load necessary libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(SoupX)
library(scDblFinder)

##### GSM7021871 Data Analysis #####
###### HCC Data Analysis ######

# Matrices files
hcc_counts <- Read10X_h5("./HCC/filtered_feature_bc_matrix.h5")
raw_hcc_counts <- Read10X_h5("./HCC/raw_feature_bc_matrix.h5")
# Create a standard Seurat Object
hcc_seurat <- CreateSeuratObject(counts = hcc_counts, project = "HCC_SingleCell")
raw_hcc_seurat <- CreateSeuratObject(counts = raw_hcc_counts, project = "HCC_SingleCell_Raw")

# Data analysis

#DimPlot(filtered_hcc_seurat, reduction = "umap", label = TRUE) + ggtitle("UMAP of Filtered HCC Data")

# Use SoupX to estimate and correct for ambient RNA contamination
soup_channel <- SoupChannel(raw_hcc_counts, hcc_counts)

# Quick clustering is required by SoupX to identify cell-specific vs ambient expression
set.seed(123)
hcc_seurat <- SCTransform(hcc_seurat) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.5) %>%
  RunUMAP(dims = 1:20)

DimPlot(hcc_seurat, reduction = "umap", label = TRUE) + ggtitle("UMAP of HCC Data for SoupX Clustering")

meta <- hcc_seurat@meta.data
soup_channel <- setClusters(soup_channel, setNames(meta$seurat_clusters, rownames(meta)))

soup_channel <- autoEstCont(
  soup_channel,
  tfidfMin = 0.1,
  soupQuantile = 0.5
)

corrected_counts <- adjustCounts(soup_channel)
# Create a new Seurat object with corrected counts
corrected_hcc_seurat <- CreateSeuratObject(counts = corrected_counts, project = "HCC_SingleCell_Corrected")

## Quality Control
corrected_hcc_seurat[["percent.mt"]] <- PercentageFeatureSet(corrected_hcc_seurat, pattern = "^MT-")
corrected_hcc_seurat[["percent.ribo"]] <- PercentageFeatureSet(corrected_hcc_seurat, pattern = "^RPS|^RPL")

# Visualize QC metrics
qc_plot <- VlnPlot(corrected_hcc_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)
print(qc_plot)

# Filter cells based on QC metrics
filtered_hcc_seurat <- subset(corrected_hcc_seurat,
  subset = nFeature_RNA > 200
  & nFeature_RNA < 2500
  & percent.mt < 5)

# Check for ambient RNA contamination
# Initial clustering to identify potential empty droplets
filtered_hcc_seurat <- NormalizeData(filtered_hcc_seurat) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:10) %>%
  FindClusters(resolution = 0.5) %>%
  RunUMAP(dims = 1:10)

DimPlot(filtered_hcc_seurat, reduction = "umap", label = TRUE) + ggtitle("UMAP of Filtered HCC Data")

# Doublet Detection using scDblFinder
sce <- scDblFinder(SingleCellExperiment(list(counts = filtered_hcc_seurat@assays$RNA$counts)))
filtered_hcc_seurat$scDblFinder.score <- sce$scDblFinder.score
filtered_hcc_seurat$scDblFinder.class <- sce$scDblFinder.class

table(filtered_hcc_seurat$scDblFinder.class)

# Visualize doublet scores
VlnPlot(filtered_hcc_seurat, group.by = "scDblFinder.class",
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')

# Remove predicted doublets
final_hcc_seurat <- subset(filtered_hcc_seurat, subset = scDblFinder.class == "singlet")
# Final clustering and visualization
final_hcc_seurat <- NormalizeData(final_hcc_seurat) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.5) %>%
  RunUMAP(dims = 1:20)

DimPlot(final_hcc_seurat, reduction = "umap", label = TRUE) + ggtitle("UMAP of Final HCC Data After QC and Doublet Removal")

