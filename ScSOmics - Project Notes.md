# Single Cell and Spatial Omics: HCC, PanIN, and PDAC Analysis

## Paper Overview and Context
### Biological Focus
The paper focuses on the transition of healthy cells into diseased states (e.g., cancer) and the spatial architecture of the **Tumor Microenvironment (TME)**. Understanding the TME is critical for identifying immune evasion and therapy resistance.

### Technical Environment
* **OS:** Artix Linux (i5-6500 + GTX 970 / i5-8265u + MX250)
* **Stack:** KDE Plasma (Wayland), Konsole, Fish shell.
* **Analysis Tools:** Seurat (R).

---

## Methodology: Single-Cell Pipeline
### 1. Single-Cell Pipeline (Cell-Centric)

Since Visium spots (55 µm) are not truly single-cell (often containing 1–10 cells), this stage focuses on treating each spot as an observation to define the transcriptomic landscape.

- **Normalization & Dimensionality Reduction:** Use `SCTransform` (Seurat) or `LogNormalize` (Scanpy). Run PCA followed by UMAP/t-SNE to visualize clusters.
    
- **Clustering:** Use the Louvain or Leiden algorithm. Your current worktree suggests you are already testing various _k_-means and graph-based clusters.
    
- **Deconvolution:** Since the spots are mixtures, use tools like **RCTD** (Robust Cell Type Decomposition), **Cell2location**, or **Tangram** to estimate the proportion of specific cell types (e.g., CAFs, epithelial cells, immune cells) within each spot using a scRNA-seq reference.
    

### 2. Spatial Analysis (Tissue-Cell Relationship)

This stage bridges the gap between expression and histology.

- **Spatially Variable Genes (SVGs):** Identify genes whose expression is not just "high" but specifically organized in space. Use **nnSVG**, **SpatialDE**, or **SPARK-X**.
    
- **Spatial Domain Detection:** Use **BayesSpace** or **SpaGCN**. These tools use both the gene expression and the $(x, y)$ coordinates to cluster spots, ensuring that "domains" (like a tumor nest or a fibrotic area) are spatially contiguous.
    
- **Cell-Cell Proximity/Interaction:** Use **Squidpy** to calculate spatial neighborhood graphs. You can perform "Permutation tests for neighborhood enrichment" to see if CAFs are significantly closer to high-grade PanIN cells than expected by chance.
    

### 3. Additional Analysis: Cell Trajectory

In the context of cancer progression (e.g., Normal $\rightarrow$ PanIN $\rightarrow$ PDAC), trajectory inference helps visualize the "pseudotime" of malignant transformation.

- **Lineage Inference:** Use **Slingshot** (R) or **PAGA** (Python). PAGA is particularly strong for spatial data as it can reconcile the discrete clusters with continuous transitions.
    
- **Spatial Trajectory:** Use **CASCAT** or **stLearn**. These tools allow you to project the pseudotime directly onto the tissue coordinates, showing exactly where in the tissue the "transition" from low-grade to high-grade PanIN is occurring.
    
- **RNA Velocity:** If you have the raw `.bam` files, use **velocyto** or **scVelo** to infer the direction of cell state changes based on spliced/unspliced mRNA ratios.

---

## Practical Applications: Multi-dataset Analysis
The current analysis involves three distinct tissue types with specific directory structures for clustering, PCA, and UMAP:

### 1. Hepatocellular Carcinoma (HCC)
* **Status:** Clustering and DiffExp performed for $k=2$ through $k=10$.
* **Methodology:** Employs both `graphclust` and `kmeans`.

### 2. Pancreatic Intraepithelial Neoplasia (PanIN)
* **Status:** Pre-processing complete (PCA/UMAP/TSNE generated).
* **Focus:** Early stage pancreatic cancer lesions.

### 3. PDAC Lymph Node (PDAClymphnode)
* **Status:** Metastatic Pancreatic Ductal Adenocarcinoma analysis.
* **Goal:** Understand spatial distribution of metastatic cells in lymphatic tissue.