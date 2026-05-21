### **Slide 1: Title Slide**
"Good morning, professors and colleagues. Today, Chiara and I are excited to present our Single Cell Spatial Omics project. Our work builds upon the foundational framework established by Atul Deshpande and colleagues in their landmark paper on mapping the spatial microenvironment landscape through latent spaces. Using high-resolution spatial transcriptomics data, our goal has been to disentangle the intricate cellular networks and physical tumor-immune interfaces that drive disease progression across different malignant architectures."

### **Slide 2: Paper Overview & Project Objectives**
"Let’s address the core methodological problem: standard 10x Genomics Visium captures spots that are 55 micrometers in diameter. Biologically, these spots do not represent single cells; they function as 'mini-bulk' transcriptomic mixtures. If you apply a global unsupervised analysis to an integrated dataset, the massive variance between organs completely squashes local biology. Therefore, the explicit objectives of our analysis were threefold: first, to isolate tissue-specific subsets and remove artificial global boundaries ; second, to apply localized variance stabilization using tailored `SCTransform` loops ; and finally, to resolve these mini-bulk spot mixtures by mapping single-cell reference atlases directly onto the slides via cellular deconvolution."

### **Slide 3: Preprocessing & Bioinformatics Pipeline**

- **Visuals:** Workflow chart showing Quality Control $\rightarrow$ SCTransform Normalization $\rightarrow$ PCA $\rightarrow$ UMAP. Divided into Single-Cell POV and Spatial POV.
"To execute this, we designed a parallel pipeline. On the data processing side, we performed a stringent data cleanup by subsetting spots based on strict mitochondrial expression thresholds. We then applied regularized negative binomial regression via `SCTransform` to stabilize variance without over-normalizing subtle cellular transitions. For dimensionality reduction, we calculated Principal Components to embed the spots into a lower-dimensional UMAP space. This allows us to examine the data from two distinct angles: the Single-Cell POV, focused on trajectory inference and spot deconvolution , and the Spatial POV, focusing on mapping these continuous cell fractions back to their physical coordinates."

### **Slide 4: Full Dataset UMAP**

- **Visuals:** Combined UMAP showing three distinct, isolated clusters colored by sample origin: HCC (Red), PanIN (Green), and PDAC (Blue).
    
- **Script:**
    

> "This is our global, merged dataset UMAP containing all captured spots across all three tissue types. Notice how the data splits into three completely isolated islands based strictly on the organ of origin. This visualization perfectly validates our core hypothesis: if we were to cluster or analyze this merged dataset globally, the macro-level transcriptomic divergence between liver tissue and pancreatic tissue would completely dominate our mathematical variance calculation. It would mask the subtle, local cellular gradients—such as localized immune infiltration—that we actually care about. This absolute separation mathematically justifies our choice to subset the data and run localized downstream pipelines."

### **Slide 5: Tissue-Specific Subset UMAPs**

- **Visuals:** Side-by-side unintegrated UMAP plots for individual subsets: UMAP of PDAC, UMAP of PanIN, and UMAP of HCC.
    
- **Script:**
    

> "By breaking the global matrix apart into tissue-specific objects, we can successfully reveal the internal sub-cluster topologies unique to each condition. Here, you can see the independent UMAP spaces for PDAC, PanIN, and HCC. Unmasking these local structures allows us to identify localized transcriptomic shifts and set up the foundation for single-cell trajectory tracking and accurate cell-type mapping."

### **Slide 6: Trajectory Analysis: HCC**

- **Visuals:** Monocle3 Trajectory Analysis plots of the HCC subset, showing cells arranged by cluster boundaries and a continuous pseudotime color gradient.
    
- **Script:**
    

> "To understand the evolutionary dynamics within these microenvironments, we ran pseudotime trajectory inference using Monocle3. Looking at the Hepatocellular Carcinoma subset, the algorithm maps a clear, continuous lineage branching out from a central trunk into distinct terminal fates. Rather than viewing these clusters as static, isolated blocks, pseudotime allows us to model the progressive transcriptomic shifts occurring as cells undergo metabolic rewiring and interact with the surrounding stroma."

### **Slide 7: Trajectory Analysis: Pancreatic Axis**

- **Visuals:** Monocle3 trajectory plot mapping the continuum between the pre-malignant PanIN and fully invasive PDAC slots.
    
- **Script:**
    

> "We applied this same trajectory framework to the pancreatic axis, linking pre-malignant PanIN modifications to invasive PDAC progression. This reveals a continuous molecular bridge rather than discrete categorical states. By tracking this axis, we can observe how the transformation of normal ductal epithelium into a pre-tumoral state dynamically correlates with changes in the expression of structural and invasive markers along the pseudotime timeline."

### **Slide 8: HCC Deconvolution Results**

- **Visuals:** Deconvolution matrix heatmap or spatial overlay showing predicted fractions of hepatocytes, immune subsets, and cholangiocarcinoma markers mapped back to the HCC slide.
    
- **Script:**
    

> "Now, let’s resolve the mini-bulk spot mixture problem using anchor-based cellular deconvolution. By projecting a single-cell liver atlas onto our HCC spatial query dataset, we successfully estimated the continuous cell-type fractions for every individual spot. The deconvolution mapping reveals a clear structural partition: we capture the distinct gradients of malignant hepatocytes, clear stromal boundaries, and localized inflammatory niches where immune cells aggregate. This bridges the gap between raw unsupervised clustering and true, highly resolved biological identity."

### **Slide 9: PanIN Deconvolution Results**

- **Visuals:** Continuous spatial probability mapping or dot plots showing normal vs. pre-tumoral ductal cell signatures distributed across the PanIN slide architecture.
    
- **Script:**
    

> "Moving to the pre-malignant pancreatic model, the deconvolution analysis allows us to look closely at the ductal cell populations. The single-cell reference mapping separates healthy, normal ductal signatures from pre-tumoral epithelial variants. Crucially, when we look at the continuous cell fractions, we can physically map where these healthy ductal cells begin losing their normal epithelial signatures and start acquiring early pre-tumoral traits within the developing lesions. This demonstrates the power of deconvolution to capture early-stage transformation gradients before clear histological changes are visible."

### **Slide 10: PDAC Deconvolution Results**

- **Visuals:** Spatial deconvolution profiles displaying malignant pancreatic ductal cell fractions completely surrounded by localized immune cell aggregates.
    
- **Script:**
    

> "Finally, we look at the fully developed, invasive PDAC sample. The deconvolution vectors here reveal a stark, highly specialized microenvironmental layout. The malignant pancreatic ductal cell signatures are heavily concentrated within dense tumor cores. However, by mapping continuous cell fractions, we show that these malignant clusters are physically surrounded by and locked into complex immune cell niches. This visualizes the classic desmoplastic and immunosuppressive boundary characteristic of pancreatic cancer, proving that our workflow can accurately map complex single-cell interactions directly onto tissue coordinates."

### **Slide 11: Conclusion & Key Takeaways**

- **Visuals:** Final summary slide highlighting: 1. Success of localized processing over global workflows. 2. Resolution of mini-bulk spots via anchor transfer. 3. Mapping of continuous disease progression axes.
    
- **Script:**
    

> "In conclusion, our analysis demonstrates that breaking away from global variance metrics is mathematically necessary to capture fine-grained tissue biology. By implementing localized `SCTransform` workflows and leveraging anchor-based single-cell deconvolution, we successfully turned mixed, mini-bulk Visium spots into interpretable, continuous landscapes of interacting cell types. We have tracked these dynamics from early pre-malignancy in PanIN up to invasive PDAC architectures. Thank you very much for your time, and Chiara and I are now open to any questions you may have."

