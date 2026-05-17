# Libraries
library(Signac)
library(tidyverse)
library(SeuratWrappers)
library(Seu)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)

##### Trajectory analysis for HCC; PanIN and PDAC tissues #####
# 1. Convert Seurat Object to Monocle3 Cell Data Set (CDS)
hcc <- CreateSeuratObject(counts = Read10X_h5("HCC/filtered_feature_bc_matrix.h5"))
panin <- CreateSeuratObject(counts = Read10X_h5("PanIN/filtered_feature_bc_matrix.h5"))
pdac <- CreateSeuratObject(counts = Read10X_h5("PDAClymphnode/filtered_feature_bc_matrix.h5"))

merged_seurat <- merge(hcc, y = c(panin, pdac), add.cell.ids = c("HCC", "PanIN", "PDAC"), project = "TrajectoryAnalysis")

rm(hcc, panin, pdac) # Clean up memory

# Quick QC and normalization (you can replace this with your existing Seurat workflow)

cleaned_seurat <- merged_seurat %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.5)
  
# Mithochndrial RNA filtering
cleaned_seurat <- cleaned_seurat %>%
  PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
  subset(subset = nCount_RNA > 200 &
         nFeature_RNA < 2500 &
         percent.mt < 5)

# This retains your normalized counts and metadata (including your SoupX cleaned counts)
cleaned_seurat <- JoinLayers(cleaned_seurat)
cds <- as.cell_data_set(cleaned_seurat, assay = "RNA")


# 2. Synchronize Reductions
# Monocle3 expects reductions to be explicitly mapped for graph learning.
# We copy the UMAP coordinates from the Seurat slot to the Monocle3 slot.
reducedDims(cds)$UMAP <- Embeddings(cleaned_seurat, reduction = "umap")

# 3. Synchronize Cluster Assignments
# Monocle3 requires its own internal clustering format before learning a graph.
# We force Monocle3 to use your Seurat clusters and assign them to a single partition.
list_cluster <- Idents(cleaned_seurat)
names(list_cluster) <- colnames(cds)

# Assign all cells to Partition 1 so Monocle tries to connect them in one graph
cds@clusters$UMAP$partitions <- factor(rep(1, ncol(cds)), levels = 1)
names(cds@clusters$UMAP$partitions) <- colnames(cds)

# Inject your Seurat clusters
cds@clusters$UMAP$clusters <- list_cluster

# 4. Learn the Graph using the imported UMAP coordinates
# We omit preprocess_cds and reduce_dimension to preserve Seurat's exact layout
cds <- learn_graph(cds, use_partition = FALSE)

# 5. Calculate Pseudotime (Interactive or Root-Based)
# Open the interactive viewer to select the root node manually on your Seurat UMAP layout
cds <- order_cells(cds)

# 6. Plot the trajectory line directly on top of your Seurat UMAP
plot_cells(cds,
           color_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_leaves = TRUE,
           label_branch_points = TRUE,
           graph_label_size = 4)

plotly::ggplotly() # Make it interactive to explore branches and pseudotime
