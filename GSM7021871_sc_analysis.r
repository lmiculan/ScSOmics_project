# Single Cell Spatial Omics Analysis for GSM7021871
library(Seurat)
library(tidyverse)
library(ggplot2)
library(Matrix)
library(data.table)
library(patchwork)
library(SeuratWrappers)
library(monocle3)
library(harmony)

# qc_pipeline(): QC/normalization + clustering + UMAP for a Seurat object.
# parameters: dimensions controls PCA/UMAP depth; res controls cluster granularity;
# assay_type selects the input assay for SCTransform; pcadim sets number of PCs computed.
qc_pipeline <- function(seurat_obj, dimensions = 20, res = 0.5, assay_type = "Spatial", pcadim = 50) {
  seurat_obj <- SCTransform(seurat_obj, assay = assay_type) |>
  RunPCA(assay = "SCT", npcs = pcadim) |>
  FindNeighbors(dims = 1:dimensions) |>
  FindClusters(resolution = res) |>
  RunUMAP(dims = 1:dimensions)
}

if (!dir.exists("results")) dir.create("results", recursive = TRUE)

set.seed(123)

hcc_data <- LoadSeuratRds("hcc_data.rds")
panin_data <- LoadSeuratRds("PanIN_data.rds")
pdac_data <- LoadSeuratRds("PDAC_data.rds")

hcc_data$orig.ident <- "HCC"
panin_data$orig.ident <- "PanIN"
pdac_data$orig.ident <- "PDAC"

# Merge samples to compare across tissues in a shared embedding.
seu_list_sub <- list(hcc_data, panin_data, pdac_data)
merged_seurat <- Reduce(function(x, y) merge(x, y), seu_list_sub)
merged_seurat <- JoinLayers(
  object = merged_seurat,
  assay = "Spatial",
  layers = "counts",
  new.layer = "counts"
)

# Set orig.ident and colours for consistent plotting across merged and individual analyses.
merged_seurat$orig.ident <- merged_seurat$orig.ident

unified_palette <- c(
  "HCC"   = "#E66101",
  "PanIN" = "#B2ABD2",
  "PDAC"  = "#5E3C99"
)

gc()
merged_seurat <- qc_pipeline(merged_seurat, dimensions = 20, res = 0.5)
DimPlot(merged_seurat, reduction = "umap", label = T, group.by = "orig.ident", cols = unified_palette, pt.size = 3) + ggtitle("UMAP of Merged Samples Colored by Sample")
ggsave(filename = "results/Merged_UMAP.png", width = 12, height = 11)
gc()

# Focus on the pancreatic disease axis by restricting to PanIN and PDAC.
pancreatic_axis <- subset(merged_seurat, subset = orig.ident %in% c("PanIN", "PDAC"))
rm(merged_seurat)

# Deconvolution: use a pancreas atlas reference to transfer cell-type labels into spatial queries.

sparse_counts <- readMM("annots/PDAC_annot/Exp_data_UMIcounts.mtx")
metadata <- fread("annots/PDAC_annot/Cells.csv", data.table = FALSE)
features <- fread("annots/PDAC_annot/Genes.txt", data.table = FALSE, header = FALSE)$V1
barcodes <- metadata$cell_name
cell_type <- metadata$cell_type

rownames(sparse_counts) <- features
colnames(sparse_counts) <- barcodes

gc()

pancreas_atlas <- CreateSeuratObject(counts = sparse_counts, meta.data = metadata, assay = "RNA")
rm(sparse_counts)

Idents(pancreas_atlas) <- metadata$cell_type

# Downsample to 100 cells to reduce compute while validating the pipeline.
set.seed(123)
pancreas_sub <- subset(pancreas_atlas, downsample = 100)

# SCTransform on the reference so anchors use SCT-normalized features.
pancreas_sub <- SCTransform(pancreas_sub, assay = "RNA", verbose = FALSE)

empty_ct <- pancreas_sub$cell_type == "" | is.na(pancreas_sub$cell_type)

sum(empty_ct)

pancreas_sub <- subset(pancreas_sub, cells = colnames(pancreas_sub)[!empty_ct])

DefaultAssay(pdac_data) <- "SCT"

# Use SCT normalization with pcaproject for stable transfer; dims=1:50 to capture broad variance.
anchors <- FindTransferAnchors(
  reference = pancreas_sub,
  query = pdac_data,
  normalization.method = "SCT",
  reference.assay = "SCT",
  query.assay = "SCT",
  reduction = "pcaproject",
  dims = 1:50
)

predictions <- TransferData(
  anchorset = anchors,
  refdata = pancreas_sub$cell_type,
  prediction.assay = TRUE,
  weight.reduction = pdac_data[["pca"]],
  dims = 1:50
)

pdac_data[["predictions"]] <- predictions
DefaultAssay(pdac_data) <- "predictions"

print(GetAssayData(pdac_data, assay = "predictions", layer = "data")[1:5, 1:5])

cell_types <- unique(rownames(pdac_data[["predictions"]]))
for (ct in cell_types) {
  p <- FeaturePlot(pdac_data, features = ct, pt.size = 3) +
    ggtitle(paste("Spatial Distribution of Predicted", ct, "in PDAC Sample"))
  ggsave(filename = paste0("results/PDAC_deconv/PDAC_Predicted_", ct, ".png"), plot = p, width = 12, height =11)
}

DefaultAssay(panin_data) <- "SCT"

anchors_pan <- FindTransferAnchors(
  reference = pancreas_sub,
  query = panin_data,
  normalization.method = "SCT",
  reference.assay = "SCT",
  query.assay = "SCT",
  reduction = "pcaproject",
  dims = 1:50
)

predictions_pan <- TransferData(
  anchorset = anchors_pan,
  refdata = pancreas_sub$cell_type,
  prediction.assay = TRUE,
  weight.reduction = panin_data[["pca"]],
  dims = 1:50
)

panin_data[["predictions"]] <- predictions_pan
DefaultAssay(panin_data) <- "predictions"

cell_types_pan <- unique(rownames(panin_data[["predictions"]]))
for (ct in cell_types_pan) {
  p <- FeaturePlot(panin_data, features = ct, pt.size = 3) +
    ggtitle(paste("Spatial Distribution of Predicted", ct, "in PanIN Sample"))
  ggsave(filename = paste0("results/PanIN_deconv/PanIN_Predicted_", ct, ".png"), plot = p, width = 12, height =11)
}

sparse_counts <- readMM("annots/HCC_annot/GSE151530_matrix.mtx.gz")
features <- fread("annots/HCC_annot/GSE151530_genes.tsv.gz", data.table = FALSE, header = FALSE)$V2
barcodes <- fread("annots/HCC_annot/GSE151530_barcodes.tsv.gz", data.table = FALSE, header = FALSE)$V1
cell_types <- fread("annots/HCC_annot/GSE151530_Info.txt", data.table = FALSE, header = T)

# Remove version suffixes to align gene symbols between atlas and spatial data.
features <- gsub("\\.\\d+$", "", features)

rownames(sparse_counts) <- features
colnames(sparse_counts) <- barcodes

# Making rownames unique to avoid issues in CreateSeuratObject
dup <- duplicated(rownames(sparse_counts))
rownames(sparse_counts) <- make.unique(rownames(sparse_counts))

gc()
hcc_atlas <- CreateSeuratObject(counts = sparse_counts, assay = "RNA")
rm(sparse_counts)
gc()

Idents(hcc_atlas) <- cell_types$Type

set.seed(123)
hcc_sub <- subset(hcc_atlas, downsample = 100)

hcc_sub <- SCTransform(hcc_sub, assay = "RNA", verbose = FALSE)

# Ensure all cells have a valid cell type for transfer; remove any with empty or NA labels.
empty_ct <- hcc_sub@active.ident == "" | is.na(hcc_sub@active.ident)

sum(empty_ct)

hcc_sub <- subset(hcc_sub, cells = colnames(hcc_sub)[!empty_ct])

DefaultAssay(hcc_data) <- "SCT"

# Use VariableFeatures(hcc_sub) and dims=1:30 to reduce noise for the HCC atlas.
anchors <- FindTransferAnchors(
  reference = hcc_sub,
  query = hcc_data,
  features = VariableFeatures(hcc_sub),
  normalization.method = "SCT",
  reference.assay = "SCT",
  query.assay = "SCT",
  reduction = "pcaproject",
  dims = 1:30
)

predictions <- TransferData(
  anchorset = anchors,
  refdata = hcc_sub@active.ident,
  prediction.assay = TRUE,
  weight.reduction = hcc_data[["pca"]],
  dims = 1:30
)

hcc_data[["predictions"]] <- predictions
DefaultAssay(hcc_data) <- "predictions"

print(GetAssayData(hcc_data, assay = "predictions", layer = "data")[1:5, 1:5])

cell_types <- unique(rownames(hcc_data[["predictions"]]))
for (ct in cell_types) {
  p <- FeaturePlot(hcc_data, features = ct, pt.size = 3) +
    ggtitle(paste("Spatial Distribution of Predicted", ct, "in HCC Sample"))
  ggsave(filename = paste0("results/HCC_deconv/HCC_Predicted_", ct, ".png"), plot = p, width = 12, height =11)
}

# Rerun QC with moderate clustering (res=0.6) for subtype separation.
pancreatic_axis <- qc_pipeline(pancreatic_axis, dimensions = 20, res = 0.6)
DimPlot(pancreatic_axis, reduction = "umap", label = TRUE, group.by = "orig.ident", cols = unified_palette, pt.size = 3) + ggtitle("UMAP of PanIN and PDAC Samples")
ggsave(filename = "results/PanIN_PDAC_UMAP.png", width = 12, height = 11)

# Trajectory analysis: reuse Seurat UMAP to order cells in monocle3.
expression_matrix <- hcc_data[["SCT"]]@counts
cell_metadata <- hcc_data@meta.data

feature_metadata <- data.frame(
  gene_short_name = rownames(expression_matrix),
  row.names = rownames(expression_matrix)
)

cds_hcc <- new_cell_data_set(
  expression_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = feature_metadata
)

reducedDims(cds_hcc)$UMAP <- Embeddings(hcc_data, reduction = "umap")

# Set a single partition to enforce a continuous trajectory.
hcc_clusters <- Idents(hcc_data)
names(hcc_clusters) <- colnames(cds_hcc)
cds_hcc@clusters$UMAP$clusters <- hcc_clusters

mock_partitions <- factor(rep(1, ncol(cds_hcc)), levels = 1)
names(mock_partitions) <- colnames(cds_hcc)
cds_hcc@clusters$UMAP$partitions <- mock_partitions

cds_hcc <- learn_graph(cds_hcc, use_partition = FALSE)
cds_hcc <- order_cells(cds_hcc, reduction_method = "UMAP")

plot_cells(cds_hcc, color_cells_by = "cluster", label_groups_by_cluster = TRUE, graph_label_size = 5, cell_size = 3) +
  ggtitle("Monocle3 Trajectory Analysis of HCC Only - Clusters")
ggsave(filename = "results/HCC_Trajectory/HCC_Trajectory_Clusters.png", width = 12, height = 11)
plot_cells(cds_hcc, color_cells_by = "pseudotime", label_groups_by_cluster = TRUE, graph_label_size = 5, cell_size = 3) +
  ggtitle("Monocle3 Trajectory Analysis of HCC Only - Pseudotime")
ggsave(filename = "results/HCC_Trajectory/HCC_Trajectory_Pseudotime.png", width = 12, height = 11)

expression_matrix_pan <- panin_data[["SCT"]]@counts
cell_metadata_pan <- panin_data@meta.data
feature_metadata_pan <- data.frame(
  gene_short_name = rownames(expression_matrix_pan),
  row.names = rownames(expression_matrix_pan)
)

cds_panin <- new_cell_data_set(
  expression_matrix_pan,
  cell_metadata = cell_metadata_pan,
  gene_metadata = feature_metadata_pan
)

reducedDims(cds_panin)$UMAP <- Embeddings(panin_data, reduction = "umap")

panin_clusters <- Idents(panin_data)
names(panin_clusters) <- colnames(cds_panin)
cds_panin@clusters$UMAP$clusters <- panin_clusters
mock_partitions_pan <- factor(rep(1, ncol(cds_panin)), levels = 1)
names(mock_partitions_pan) <- colnames(cds_panin)
cds_panin@clusters$UMAP$partitions <- mock_partitions_pan

cds_panin <- learn_graph(cds_panin, use_partition = FALSE)
cds_panin <- order_cells(cds_panin, reduction_method = "UMAP")
plot_cells(cds_panin, color_cells_by = "cluster", label_groups_by_cluster = TRUE, graph_label_size = 3, cell_size = 2) +
  ggtitle("Monocle3 Trajectory Analysis of PanIN Only - Clusters")
ggsave(filename = "results/PanIN_Trajectory/PanIN_Trajectory_Clusters.png", width = 12, height = 11)
plot_cells(cds_panin, color_cells_by = "pseudotime", label_groups_by_cluster = TRUE, graph_label_size = 3, cell_size = 2) +
  ggtitle("Monocle3 Trajectory Analysis of PanIN Only - Pseudotime")
ggsave(filename = "results/PanIN_Trajectory/PanIN_Trajectory_Pseudotime.png", width = 12, height = 11)

expression_matrix_pdac <- pdac_data[["SCT"]]@counts
cell_metadata_pdac <- pdac_data@meta.data
feature_metadata_pdac <- data.frame(
  gene_short_name = rownames(expression_matrix_pdac),
  row.names = rownames(expression_matrix_pdac)
)

cds_pdac <- new_cell_data_set(
  expression_matrix_pdac,
  cell_metadata = cell_metadata_pdac,
  gene_metadata = feature_metadata_pdac
)

reducedDims(cds_pdac)$UMAP <- Embeddings(pdac_data, reduction = "umap")

pdac_clusters <- Idents(pdac_data)
names(pdac_clusters) <- colnames(cds_pdac)
cds_pdac@clusters$UMAP$clusters <- pdac_clusters
mock_partitions_pdac <- factor(rep(1, ncol(cds_pdac)), levels = 1)
names(mock_partitions_pdac) <- colnames(cds_pdac)
cds_pdac@clusters$UMAP$partitions <- mock_partitions_pdac

cds_pdac <- learn_graph(cds_pdac, use_partition = FALSE)
cds_pdac <- order_cells(cds_pdac, reduction_method = "UMAP")
plot_cells(cds_pdac, color_cells_by = "cluster", label_groups_by_cluster = TRUE, graph_label_size = 3, cell_size = 3) +
  ggtitle("Monocle3 Trajectory Analysis of PDAC Only - Clusters")
ggsave(filename = "results/PDAC_Trajectory/PDAC_Trajectory_Clusters.png", width = 12, height = 11)
plot_cells(cds_pdac, color_cells_by = "pseudotime", label_groups_by_cluster = TRUE, graph_label_size = 3, cell_size = 3) +
  ggtitle("Monocle3 Trajectory Analysis of PDAC Only - Pseudotime")
ggsave(filename = "results/PDAC_Trajectory/PDAC_Trajectory_Pseudotime.png", width = 12, height = 11)
