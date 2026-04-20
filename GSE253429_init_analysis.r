library(rhdf5)
library(Matrix)
library(Seurat)

h5_data <- h5read("GSM8020617_PN1_filtered_feature_bc_matrix.h5", "/matrix")

sparse_mtx <- sparseMatrix(
  i = h5_data$indices + 1,
  p = h5_data$indptr,
  x = as.numeric(h5_data$data),
  dims = h5_data$shape
)

raw_genes <- as.character(h5_data$features$name)
rownames(sparse_mtx) <- make.unique(raw_genes)

raw_barcodes <- as.character(h5_data$barcodes)

if (length(raw_barcodes) == ncol(sparse_mtx)) {
  colnames(sparse_mtx) <- raw_barcodes
} else {
  colnames(sparse_mtx) <- paste0("Cell_", 1:ncol(sparse_mtx))
  warning("Dimensioni barcodes non corrispondenti. Usati nomi generici.")
}

seurat_obj <- CreateSeuratObject(counts = sparse_mtx, project = "Pancreas")

##### Analysis ######

saveRDS(seurat_obj, file = "seurat_pancreas.rds")

# QC
print(seurat_obj)

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

violin_plot <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 2)
violin_plot

# Subsetting
seurat_obj_sub <- subset(seurat_obj, subset = nFeature_RNA > 200 &
                       nCount_RNA > 2500 &
                       percent.mt < 5)
print(seurat_obj_sub)


# Normalization
seurat_obj_sub <- SCTransform(seurat_obj_sub, vars.to.regress = "percent.mt", verbose = FALSE)

# Dimension reduction
seurat_obj_sub <- RunPCA(seurat_obj_sub, verbose = FALSE)
seurat_obj_sub <- RunUMAP(seurat_obj_sub, dims = 1:30, verbose = F)
seurat_obj_sub <- RunTSNE(seurat_obj_sub, dims = 1:30, verbose = FALSE)

DimPlot(seurat_obj_sub, reduction = "umap", label = TRUE) + NoLegend()
DimPlot(seurat_obj_sub, reduction = "tsne", label = TRUE) + NoLegend()

# Clustering
seurat_obj_sub <- FindNeighbors(seurat_obj_sub, dims = 1:30, verbose = FALSE)
seurat_obj_sub <- FindClusters(seurat_obj_sub, resolution = 0.5, verbose = FALSE)

DimPlot(seurat_obj_sub, reduction = "umap", label = TRUE) + NoLegend()

# Marker genes
markers <- FindAllMarkers(seurat_obj_sub, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
top2 <- markers %>% group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)
DoHeatmap(seurat_obj_sub, features = top2$gene) + NoLegend()
