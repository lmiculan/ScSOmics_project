library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)

##### HCC Data Analysis #####
# Da cambiare
hcc_data <- Load10X_Spatial(data.dir = "HCC/")

#2. Add og clusters to metadata
original_labels <- read.delim("HCC/analysis/clustering/graphclust/clusters.csv", sep = ",")
original_labels <- original_labels %>%
  select(Barcode, Cluster) %>%
  rename(orig.ident = Barcode, original_clusters = Cluster)

hcc_data <- AddMetaData(hcc_data, metadata = original_labels)

##### Quality control #####

# Initial processing
hcc_data <- SCTransform(hcc_data, assay = "Spatial", verbose = FALSE)
hcc_data <- RunPCA(hcc_data, verbose = FALSE)
hcc_data <- FindNeighbors(hcc_data, dims = 1:30)
hcc_data <- FindClusters(hcc_data, resolution = 0.3)

# Visualization for comparison
p1 <- SpatialDimPlot(hcc_data, group.by = "seurat_clusters")
p2 <- SpatialDimPlot(hcc_data, group.by = "original_clusters")
p1 + p2

library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)

##### PanIN Data Analysis #####

panin_data <- Load10X_Spatial(data.dir = "PanIN/")

#2. Add og clusters to metadata
original_labels <- read.delim("PanIN/analysis/clustering/1152735_analysis__clustering_kmeans_10_clusters_clusters.csv", sep = ",")
original_labels <- original_labels %>%
  select(Barcode, Cluster) %>%
  rename(orig.ident = Barcode, original_clusters = Cluster)

panin_data <- AddMetaData(panin_data, metadata = original_labels)

##### Quality control #####

# Initial processing
panin_data <- SCTransform(panin_data, assay = "Spatial", verbose = FALSE)
panin_data <- RunPCA(panin_data, verbose = FALSE)
panin_data <- FindNeighbors(panin_data, dims = 1:30)
panin_data <- FindClusters(panin_data, resolution = 0.3)

# Visualization for comparison
p1 <- SpatialDimPlot(panin_data, group.by = "seurat_clusters")
p2 <- SpatialDimPlot(panin_data, group.by = "original_clusters")
p1 + p2
