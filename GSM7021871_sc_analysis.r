# Single Cell Spatial Omics Analysis for GSM7021871
# Load necessary libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(SoupX)
library(scDblFinder)
library(SeuratWrappers)
library(monocle3)
library(dplyr)

##### Functions #####
# Quality control pipeline for Seurat objects
qc_pipeline <- function(seurat_obj, dimensions, res) {
  seurat_obj <- SCTransform(seurat_obj) |>
  FindVariableFeatures() |>
  ScaleData() |>
  RunPCA() |>
  FindNeighbors(dims = 1:dimensions) |>
  FindClusters(resolution = res) |>
  RunUMAP(dims = 1:dimensions)
}

# SoupX function on multi sample Seurat
run_soupx_automated <- function(seurat_obj, raw_mtx, filtered_mtx) {
  print(paste("Running SoupX on", seurat_obj@project.name))
  
  # 1. Access Gene Names specifically from the Dimnames list
  raw_genes <- raw_mtx@Dimnames[[1]]
  filt_genes <- filtered_mtx@Dimnames[[1]]
  
  print(paste("Raw matrix gene dimensions:", length(raw_genes)))
  print(paste("Filtered matrix gene dimensions:", length(filt_genes)))
  
  # 2. Intersect ONLY the gene names
  common_genes <- intersect(raw_genes, filt_genes)
  
  if(length(common_genes) == 0) {
     stop("Zero common genes found. Check if one matrix uses Ensembl IDs and the other Symbols.")
  }

  # 3. Subset matrices using the character vector of gene names
  raw_aligned <- raw_mtx[common_genes, ]
  filt_aligned <- filtered_mtx[common_genes, ]

  # 4. Create SoupChannel
  sc <- SoupChannel(raw_aligned, filt_aligned)
  
  # 5. Handle Barcode Mismatches
  clusters <- Idents(seurat_obj)
  names(clusters) <- gsub("-1_.*", "-1", names(clusters))
  
  # Ensure barcodes match the aligned filtered matrix
  # filt_aligned@Dimnames[[2]] contains the barcodes
  common_barcodes <- intersect(names(clusters), colnames(filt_aligned))
  clusters <- clusters[common_barcodes]
  
  # 6. Run SoupX Pipeline
  sc <- setClusters(sc, clusters)
  sc <- autoEstCont(sc, doPlot = FALSE, forceAccept=TRUE)
  out <- adjustCounts(sc)
  
  # 7. Create cleaned Seurat Object
  # Subsetting the original metadata to match surviving barcodes
  cleaned_meta <- seurat_obj@meta.data[common_barcodes, ]
  cleaned_seurat <- CreateSeuratObject(counts = out, 
                                       meta.data = cleaned_meta, 
                                       project = seurat_obj@project.name)

  # Clean up memory immediately
  rm(raw_mtx, filtered_mtx, raw_aligned, filt_aligned, out); gc()

  return(cleaned_seurat)
}

# Multiple sample doublet detection using scDblFinder
dblfinder <- function(seurat_obj) {
  counts_mat <- LayerData(seurat_obj[["RNA"]], layer = "counts")

  sce <- SingleCellExperiment(
    assays = list(counts = counts_mat),
    colData = seurat_obj@meta.data
  )

  sce <- scDblFinder(sce, samples = "orig.ident")

  seurat_obj$scDblFinder.score <- sce$scDblFinder.score
  seurat_obj$scDblFinder.class <- sce$scDblFinder.class
  seurat_obj
}


##### GSM7021871 Data Analysis #####

# Directory for results
if (!dir.exists("results")) dir.create("results", recursive = TRUE)

###### Data Analysis ######

filtered_data_dirs <- list("./HCC/filtered_feature_bc_matrix.h5", "./PanIN/filtered_feature_bc_matrix.h5", "./PDAClymphnode/filtered_feature_bc_matrix.h5")
raw_data_dirs <- list("./HCC/raw_feature_bc_matrix.h5", "./PanIN/raw_feature_bc_matrix.h5", "./PDAClymphnode/raw_feature_bc_matrix.h5")

# Matrices files
hcc_counts <- Read10X_h5(filtered_data_dirs[[1]])
raw_hcc_counts <- Read10X_h5(raw_data_dirs[[1]])

panin <- Read10X_h5(filtered_data_dirs[[2]])
raw_panin <- Read10X_h5(raw_data_dirs[[2]])

pdac <- Read10X_h5(filtered_data_dirs[[3]])
raw_pdac <- Read10X_h5(raw_data_dirs[[3]])

filtered_list <- list(hcc_counts, panin, pdac)
raw_list <- list(raw_hcc_counts, raw_panin, raw_pdac)

# Create a standard Seurat Object
hcc_seurat <- CreateSeuratObject(counts = hcc_counts, project = "HCC_SingleCell")
raw_hcc_seurat <- CreateSeuratObject(counts = raw_hcc_counts, project = "HCC_SingleCell_Raw")

panin_seurat <- CreateSeuratObject(counts = panin, project = "PanIN_SingleCell")
raw_panin_seurat <- CreateSeuratObject(counts = raw_panin, project = "PanIN_SingleCell_Raw")

pdac_seurat <- CreateSeuratObject(counts = pdac, project = "PDAC_SingleCell")
raw_pdac_seurat <- CreateSeuratObject(counts = raw_pdac, project = "PDAC_SingleCell_Raw")

seu_list <- list(hcc_seurat, panin_seurat, pdac_seurat)

raw_seu_list <- list(raw_hcc_seurat, raw_panin_seurat, raw_pdac_seurat)

# Initial QC for SoupX
seu_list_qc <- map(seu_list, ~ qc_pipeline(.x, dimensions = 20, res = 0.5))

soupx_list <- pmap(
  list(seu_list_qc, raw_list, filtered_list),
  ~ run_soupx_automated(..1, ..2, ..3)
)
# Second quality control
soupx_list <- map(soupx_list, ~ {
  .x[["percent.mt"]] <- PercentageFeatureSet(.x, pattern = "^MT-")
  .x[["percent.ribo"]] <- PercentageFeatureSet(.x, pattern = "^RPS|^RPL")
  return(.x)
})

# Plot SoupX results for each sample
lapply(seq_along(soupx_list), function(i) {
  p <- VlnPlot(soupx_list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) +
    ggtitle(paste("Percentage of Mitochondrial DNA for", soupx_list[[i]]@project.name))
  print(p)
})

# Subset based on QC metrics
soupx_list_sub <- map(soupx_list, ~ subset(.x, subset =
                                            nFeature_RNA > 200 &
                                            nFeature_RNA < 2500 &
                                            percent.mt < 5))

# Merging all samples into a single Seurat object for downstream analysis
merged_seurat <- Reduce(function(x, y) merge(x, y), soupx_list_sub)
merged_seurat <- JoinLayers(
  object = merged_seurat,
  assay = "RNA",
  layers = "counts",
  new.layer = "counts"
)

# Delete single dataset to retrive memory
rm(list = c("hcc_seurat", "panin_seurat", "pdac_seurat", "raw_panin_seurat", "raw_hcc_counts", "raw_pdac_seurat"))

##### Doublet Detection using scDblFinder #####
merged_seurat <- dblfinder(merged_seurat)

# Subsetting to keep only singlets
merged_seurat <- subset(merged_seurat, subset = scDblFinder.class == "singlet")

# Final clustering and visualization
final_seurat <- qc_pipeline(merged_seurat, dimensions = 20, res = 0.5)
ElbowPlot(final_seurat, ndims = 50) + ggtitle("Elbow Plot for PCA - Final Seurat Object")

# Dimplots for visualiztion
DimPlot(final_seurat, reduction = "umap", group.by = "orig.ident",label = TRUE) + ggtitle("UMAP of Final Data After QC and Doublet Removal - Samples")

##### Differentially expressed analysis #####

tissues <- c("HCC_SingleCell", "PanIN_SingleCell", "PDAC_SingleCell")

# QC once per tissue
seu_qc <- tissues |>
  setNames(tissues) |>
  lapply(function(t) {
    seu_t <- subset(final_seurat, subset = orig.ident == t)
    qc_pipeline(seu_t, dimensions = 20, res = 0.5)
  })

run_de_on_qc <- function(seu_t, tissue) {
  DefaultAssay(seu_t) <- "SCT"
  if (is.null(Idents(seu_t))) Idents(seu_t) <- seu_t$seurat_clusters

  markers <- FindAllMarkers(
    seu_t,
    assay = "SCT",
    slot = "data",
    test.use = "wilcox",
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )

  top5 <- markers |>
    slice_max(avg_log2FC, n = 5, with_ties = FALSE, by = clusters)

  top_genes <- unique(top5$gene)

  heatmap <- DoHeatmap(
    seu_t,
    features = intersect(top_genes, rownames(seu_t)),
    assay = "SCT",
    slot = "scale.data"
  ) + NoLegend() + ggtitle(paste("Top 5 markers:", tissue))

  feature_plot <- if ("Spatial" %in% names(seu_t@assays)) {
    SpatialFeaturePlot(seu_t, features = top_genes, ncol = 3)
  } else {
    FeaturePlot(seu_t, features = top_genes, ncol = 3)
  }

  list(markers = markers, top5 = top5, heatmap = heatmap, feature_plot = feature_plot)
}

de_results <- Map(run_de_on_qc, seu_qc, names(seu_qc))
de_results
 
## Trajectory analysis with pseudotime ordering ##
# 7. Order cells in pseudotime (root cells can be specified if known, otherwise Monocle3 will attempt to infer them)
cds_pancreatic <- order_cells(cds_pancreatic)
plot_cells(cds_pancreatic, color_cells_by = "pseudotime", label_groups_by_cluster = TRUE, graph_label_size = 3, cell_size = 1) +
  ggtitle("Monocle3 Trajectory Analysis with Pseudotime Coloring")

# 1. Isolate only the HCC cells
hcc_only <- subset(final_seurat, subset = orig.ident %in% "HCC_SingleCell")

# 2. Extract the matrix and metadata
expression_matrix <- hcc_only[["SCT"]]@counts
cell_metadata <- hcc_only@meta.data

feature_metadata <- data.frame(
  gene_short_name = rownames(expression_matrix),
  row.names = rownames(expression_matrix)
)

# 3. Build the Monocle3 object
cds_hcc <- new_cell_data_set(
  expression_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = feature_metadata
)

# 4. Inject the existing HCC UMAP coordinates
reducedDims(cds_hcc)$UMAP <- Embeddings(hcc_only, reduction = "umap")

# 5. Assign a single mock partition to ensure a continuous path
hcc_clusters <- Idents(hcc_only)
names(hcc_clusters) <- colnames(cds_hcc)
cds_hcc@clusters$UMAP$clusters <- hcc_clusters

mock_partitions <- factor(rep(1, ncol(cds_hcc)), levels = 1)
names(mock_partitions) <- colnames(cds_hcc)
cds_hcc@clusters$UMAP$partitions <- mock_partitions

# 6. Learn the specific trajectory graph for HCC
cds_hcc <- learn_graph(cds_hcc, use_partition = FALSE)
cds_hcc <- order_cells(cds_hcc, reduction_method = "UMAP")
plot_cells(cds_hcc, color_cells_by = "cluster", label_groups_by_cluster = TRUE, graph_label_size = 3, cell_size = 1) +
  ggtitle("Monocle3 Trajectory Analysis of HCC Only")

##### Archive #####

## Ensure SCT is set as the default assay if present
#DefaultAssay(final_seurat) <- "SCT"
#
## Ensure identities are set to Seurat clusters
#if (is.null(Idents(final_seurat))) Idents(final_seurat) <- final_seurat$seurat_clusters
#
## Find markers for all clusters (only positive markers)
#markers_all <- FindAllMarkers(final_seurat,
#                              assay = DefaultAssay(final_seurat),
#                              slot = "data",
#                              test.use = "wilcox",
#                              only.pos = TRUE,
#                              min.pct = 0.25,
#                              logfc.threshold = 0.25)
#
## Order and save marker table
#markers_all <- markers_all %>% arrange(orig.ident, desc(avg_log2FC))
#write.csv(markers_all, file = file.path("results", "markers_by_cluster.csv"), row.names = FALSE)
#saveRDS(markers_all, file = file.path("results", "markers_all_clusters.rds"))
#
## Select top 10 markers per cluster and plot heatmap
#top10 <- markers_all %>% group_by(orig.ident) %>% filter(p_val_adj < 0.05)%>% slice_max(order_by = avg_log2FC, n = 10) %>% ungroup()
#top_genes <- unique(top10$gene)
#
## Draw heatmap of top markers (uses scaled data from SCT if available)
#heatmap_plot <- DoHeatmap(final_seurat, features = top_genes, assay = DefaultAssay(final_seurat), slot = "scale.data") +
#  NoLegend() +
#  ggtitle("Top 10 markers per cluster")
#print(heatmap_plot)
#
### Feature scatter plots of genes ##
#### Pancreatic axis genes ###
#pancreatic_axis <- subset(final_seurat, subset = orig.ident %in% c("PanIN_SingleCell", "PDAC_SingleCell"))
#
## Rerun qc pipeline on the subset to ensure proper dimensionality reduction and clustering
#pancreatic_axis <- qc_pipeline(pancreatic_axis, dimensions = 20, res = 0.5)
#
#markers_all <- FindAllMarkers(pancreatic_axis,
#                              assay = DefaultAssay(final_seurat),
#                              slot = "data",
#                              test.use = "wilcox",
#                              only.pos = TRUE,
#                              min.pct = 0.25,
#                              logfc.threshold = 0.25)
#
#markers_all <- markers_all %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 5)
#DoHeatmap(pancreatic_axis, features = markers_all$gene[1:20], assay = DefaultAssay(pancreatic_axis), slot = "scale.data") +
#  NoLegend() +
#  ggtitle("Top 20 markers for PanIN and PDAC")
#feature_plots
### TBA cluster identification ##
#
### Trajectory analysis ##
## 1. Isolate PanIN and PDAC samples for trajectory analysis
#
## 2. Extract the SCT Assay counts from the subset
#expression_matrix <- pancreatic_axis[["SCT"]]@counts
#cell_metadata <- pancreatic_axis@meta.data
#
#feature_metadata <- data.frame(
#  gene_short_name = rownames(expression_matrix),
#  row.names = rownames(expression_matrix)
#)
#
## 3. Create the new Monocle3 object
#cds_pancreatic <- new_cell_data_set(
#  expression_matrix,
#  cell_metadata = cell_metadata,
#  gene_metadata = feature_metadata
#)
#
## 4. Inject the subsetted UMAP coordinates
## Seurat automatically subsets the embeddings when you subset the object
#reducedDims(cds_pancreatic)$UMAP <- Embeddings(pancreatic_axis, reduction = "umap")
#
## 5. Synchronize clusters for graph learning
#pan_clusters <- Idents(pancreatic_axis)
#names(pan_clusters) <- colnames(cds_pancreatic)
#cds_pancreatic@clusters$UMAP$clusters <- pan_clusters
#
## Mock partition to ensure a single unified graph
#mock_partitions <- factor(rep(1, ncol(cds_pancreatic)), levels = 1)
#names(mock_partitions) <- colnames(cds_pancreatic)
#cds_pancreatic@clusters$UMAP$partitions <- mock_partitions
#
## 6. Learn the graph on the pancreatic subset
#cds_pancreatic <- learn_graph(cds_pancreatic, use_partition = FALSE)
#
#plot_cells(cds_pancreatic, color_cells_by = "partition", label_groups_by_cluster = TRUE, graph_label_size = 3) +
#  ggtitle("Monocle3 Trajectory Analysis of PanIN and PDAC")
#