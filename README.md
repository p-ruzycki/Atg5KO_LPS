# SingleCell_RNA_Seq

## Single-Cell RNA-seq Analysis Plots

Click on the links below to view the interactive plots:

- [Plot 1](https://sohinireddy45.github.io/SingleCell_RNA_Seq/plot1.html)
- [Plot 2](https://sohinireddy45.github.io/SingleCell_RNA_Seq/plot2.html)
- [Plot 3](https://sohinireddy45.github.io/SingleCell_RNA_Seq/plot3.html)
- [Plot 4](https://sohinireddy45.github.io/SingleCell_RNA_Seq/plot4.html)

## 📊 Single-Cell RNA-seq Analysis

🔗 **View Interactive Plots:** [Click here](https://sohinireddy45.github.io/SingleCell_RNA_Seq/)

## Single-Cell RNA-seq Analysis Plots

Below are the interactive plots:

## UMAP Plot for Single-Cell RNA-seq

![UMAP Plot](https://github.com/sohinireddy45/SingleCell_RNA_Seq/blob/main/Umpa_plot_1.png?raw=true)

👉 [Click here for the full interactive version](https://sohinireddy45.github.io/SingleCell_RNA_Seq/)

### **Overview**  
This repository contains a computational pipeline for **Single-Cell RNA Sequencing (scRNA-seq) analysis**, leveraging **Seurat** and other bioinformatics tools for **quality control, clustering, differential expression analysis, and pathway enrichment**.  

### **Pipeline Summary**  
The workflow consists of the following key steps:  
1. **Preprocessing & Quality Control** – Data import, Seurat object creation, QC filtering, and normalization.  
2. **Dimensionality Reduction & Clustering** – PCA, UMAP visualization, and clustering.  
3. **Differential Expression Analysis (DEA)** – Identifying differentially expressed genes (DEGs).  
4. **Visualization** – Generating **UMAP plots, volcano plots, and heatmaps**.  
5. **Pathway Enrichment Analysis** – Using GSEA to identify biological pathways.  

---
### **1. Preprocessing & Quality Control**  

This step involves:  
- **Importing raw count matrices** from 10x Genomics  
- **Creating a Seurat object** for downstream analysis  
- **Filtering low-quality cells** based on mitochondrial content and gene count thresholds  
- **Normalizing the dataset** using the LogNormalize method  

The goal is to **remove low-quality cells and batch effects** while preserving biologically relevant variation.

#### **Code: Import Data & Create Seurat Object**  
```r
library(Seurat)

# Load count matrices
data <- Read10X(data.dir = "data/")  # Update with the correct file path

# Create a Seurat object
seurat_obj <- CreateSeuratObject(counts = data, project = "SingleCellRNAseq")

# View basic summary
print(seurat_obj)
```

---

#### **Code: Quality Control Filtering**  
Cells are filtered based on the **percentage of mitochondrial genes** and **minimum gene count** to remove low-quality cells.  

```r
# Calculate mitochondrial gene percentage
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")

# Filter cells based on mitochondrial content (<25%) and number of detected genes (>500)
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 500 & percent.mt < 25)

# Display the number of remaining cells
print(paste("Remaining cells after filtering:", ncol(seurat_obj)))
```

---

#### **Code: Normalization**  
Data is **normalized** using the **LogNormalize method**, scaling gene expression counts for comparability.  

```r
# Normalize data
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# View normalized expression levels
head(seurat_obj@assays$RNA@data)
```

---

## **2. Dimensionality Reduction & Clustering**  

This step involves:  
- **Identifying highly variable genes**  
- **Performing Principal Component Analysis (PCA)** for dimensionality reduction  
- **Determining optimal principal components (PCs) for clustering**  
- **Clustering cells** using shared nearest neighbor (SNN) modularity optimization  
- **Visualizing clusters** using Uniform Manifold Approximation and Projection (UMAP)  

---

### **Code: Identify Highly Variable Genes**  
To focus on biologically significant genes, we identify the top **2,000 most variable genes** across cells.  

```r
# Identify highly variable genes
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# View top 10 variable genes
head(VariableFeatures(seurat_obj), 10)
```

---

### **Code: Scaling & Principal Component Analysis (PCA)**  
We **scale the data** and perform **PCA** to reduce dimensionality while preserving variance.  

```r
# Scale the data
seurat_obj <- ScaleData(seurat_obj)

# Perform PCA
seurat_obj <- RunPCA(seurat_obj)

# Visualize variance explained by PCs
ElbowPlot(seurat_obj)
```
> **Interpretation**: The **Elbow Plot** helps determine the number of PCs to retain for clustering.

---

### **Code: Clustering Cells**  
We **find neighbors and cluster cells** based on the retained PCs.  
```r
# Determine optimal number of PCs (e.g., 1 to 11 based on ElbowPlot)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:11)

# Cluster cells at multiple resolutions; we use resolution = 1.0 for final clustering
seurat_obj <- FindClusters(seurat_obj, resolution = 1.0)

# View cluster assignments
table(Idents(seurat_obj))
```

---

### **Code: UMAP Visualization**  
We apply **UMAP** for **visualizing clusters** in low-dimensional space.  

```r
# Run UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:11)

# Plot UMAP with cluster labels
DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5)
```
> **Interpretation**: The **UMAP plot** shows clusters of transcriptionally distinct cell populations.

---

## **3. Differential Expression Analysis (DEA)**  

This step involves:  
- **Comparing gene expression between WT and KO conditions**  
- **Identifying differentially expressed genes (DEGs)** using the `FindMarkers()` function  
- **Performing intra-group heterogeneity analysis** (Cluster 2 vs Cluster 1)  
- **Filtering significant genes** (p_adj < 0.05, |log2FC| > 1)  
- **Saving the results for downstream analysis**  

---

### **Code: DEA - WT vs. KO in Müller Glia (Cluster 1)**  
We compare gene expression between **wild-type (WT) and knockout (KO) conditions** within **Müller glia (Cluster 1)**.

```r
# Perform DEA between WT and KO within Müller Glia (Cluster 1)
deg_results_WT_vs_KO <- FindMarkers(seurat_obj, ident.1 = "WT", ident.2 = "KO", min.pct = 0.25)

# Filter significantly differentially expressed genes (DEGs)
deg_results_WT_vs_KO <- deg_results_WT_vs_KO[deg_results_WT_vs_KO$p_val_adj < 0.05 & abs(deg_results_WT_vs_KO$log2FC) > 1, ]

# Save results
write.csv(deg_results_WT_vs_KO, "results/DEA_WT_vs_KO.csv")

# View top DEGs
head(deg_results_WT_vs_KO)
```
> **Interpretation**: The output shows **significantly upregulated and downregulated genes** in KO vs. WT.

---

### **Code: DEA - Müller Glia Subclusters (Cluster 2 vs. Cluster 1)**  
To assess **intra-group heterogeneity**, we compare expression between **Cluster 2 (disease) and Cluster 1 (control).**  

```r
# Perform DEA between Müller Glia subclusters (Cluster 2 vs Cluster 1)
deg_results_Cluster2_vs_Cluster1 <- FindMarkers(seurat_obj, ident.1 = "Cluster2", ident.2 = "Cluster1", min.pct = 0.25)

# Filter significant DEGs
deg_results_Cluster2_vs_Cluster1 <- deg_results_Cluster2_vs_Cluster1[deg_results_Cluster2_vs_Cluster1$p_val_adj < 0.05 & abs(deg_results_Cluster2_vs_Cluster1$log2FC) > 1, ]

# Save results
write.csv(deg_results_Cluster2_vs_Cluster1, "results/DEA_Cluster2_vs_Cluster1.csv")

# View top DEGs
head(deg_results_Cluster2_vs_Cluster1)
```
> **Interpretation**: This highlights gene expression differences between **Müller glia subclusters**, capturing **heterogeneity within the population**.

---

### **Next Step: Data Visualization**  
Now that we have identified DEGs, the next step is to **generate visualizations** such as:  
- **Volcano plots for DEGs**  
- **Feature plots for gene expression patterns**  
- **UMAP overlays for marker genes**  

## **4. Data Visualization**  

This step involves:  
- **Generating volcano plots** to visualize differentially expressed genes (DEGs)  
- **Creating UMAP feature plots** for key marker genes  
- **Plotting expression patterns** using dot plots and violin plots  

### **Code: Volcano Plot - WT vs. KO in Müller Glia**  
A **volcano plot** visualizes upregulated and downregulated genes between WT and KO conditions.  

```r
library(ggplot2)
library(ggrepel)

# Load DEA results
deg_results_WT_vs_KO <- read.csv("results/DEA_WT_vs_KO.csv", row.names = 1)

# Define significance thresholds
p_val_adj <- 0.05
FCcutoff <- 1

# Convert p-values to -log10 scale
deg_results_WT_vs_KO$neg_log10_pval <- -log10(deg_results_WT_vs_KO$p_val_adj)

# Categorize genes: Significant (Red) or Non-significant (Gray)
deg_results_WT_vs_KO$category <- "NS"
deg_results_WT_vs_KO$category[deg_results_WT_vs_KO$p_val_adj < p_val_adj & abs(deg_results_WT_vs_KO$log2FC) > FCcutoff] <- "Significant"

# Select key genes for labeling
selected_genes <- c("Lcn2", "Ccl2", "Lif", "Ccl7", "Gfap", "Bdnf", "Tnfaip2", "Aqp4", "Crb1", "Glul", "Rlbp1", "Kcnj10", "Car14")
selected_genes_df <- deg_results_WT_vs_KO[rownames(deg_results_WT_vs_KO) %in% selected_genes, ]

# Create Volcano Plot
p <- ggplot(deg_results_WT_vs_KO, aes(x = log2FC, y = neg_log10_pval)) +
  
  geom_point(aes(color = category), alpha = 0.75, size = 4) +  
  scale_color_manual(values = c("NS" = "gray", "Significant" = "red")) +  

  geom_point(data = selected_genes_df, aes(x = log2FC, y = neg_log10_pval), color = "black", size = 5) +

  geom_text_repel(data = selected_genes_df, aes(x = log2FC, y = neg_log10_pval, label = rownames(selected_genes_df)), 
                  size = 5, fontface = "italic", box.padding = 1, point.padding = 0.3) +

  geom_hline(yintercept = -log10(p_val_adj), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-FCcutoff, FCcutoff), linetype = "dashed", color = "black") +

  labs(title = "Volcano Plot: WT vs. KO in Müller Glia", 
       x = expression(~Log[2]~ "Fold Change"), 
       y = expression(-Log[10]~ "FDR"),
       color = "Expression") +
  
  theme_minimal()

# Save and display
ggsave("figures/volcano_WT_vs_KO.pdf", plot = p, width = 10, height = 8, units = "in", dpi = 300)
print(p)
```
> **Interpretation**:  
- **Red dots** represent significantly differentially expressed genes.  
- **Black labels** highlight key genes of interest.  
- **Dashed lines** indicate statistical significance thresholds.  

---

### **Code: UMAP Feature Plot of Selected Genes**  
We visualize the expression patterns of **key marker genes** across **Müller glia states**.

```r
# Define genes to visualize
marker_genes <- c("Crb1", "Glul", "Rlbp1", "Gfap", "Lif", "Ccl2")

# Generate feature plots
FeaturePlot(seurat_obj, features = marker_genes, cols = c("lightgray", "red"), min.cutoff = "q10", max.cutoff = "q90", ncol = 3)
```
> **Interpretation**:  
- **Dark red areas** indicate high expression of a gene.  
- **Lighter gray areas** show lower expression.  

---

### **Code: Dot Plot of Cell Type Marker Genes**  
Dot plots visualize **expression levels across cell types**.

```r
# Define cell type markers
cell_type_markers <- c("Aqp4", "Crb1", "Glul", "Rlbp1", "Gfap", "Lcn2", "Cxcl2", "Ccl2")

# Generate dot plot
DotPlot(seurat_obj, features = cell_type_markers, cols = c("blue", "red"), dot.scale = 5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```
> **Interpretation**:  
- **Dot size** represents the percentage of cells expressing a gene.  
- **Color intensity** shows expression levels.  

---

## **5. Pathway Enrichment Analysis (GSEA)**  

This step involves:  
- **Performing Gene Set Enrichment Analysis (GSEA)** to identify biological pathways affected by differentially expressed genes (DEGs).  
- **Filtering enriched pathways** based on statistical significance (NES ≥ 2, FDR q-value ≤ 0.05).  
- **Visualizing pathway enrichment results** using dot plots.  

---

### **Code: Preparing Data for GSEA**  
GSEA requires a **ranked gene list** based on log2 fold change.

```r
library(clusterProfiler)
library(org.Mm.eg.db)  # Load appropriate organism database (e.g., org.Hs.eg.db for human)
library(ggplot2)

# Load DEA results
deg_results_WT_vs_KO <- read.csv("results/DEA_WT_vs_KO.csv", row.names = 1)

# Create ranked gene list for GSEA
gene_list <- deg_results_WT_vs_KO$log2FC
names(gene_list) <- rownames(deg_results_WT_vs_KO)

# Sort genes by log2 fold change
gene_list <- sort(gene_list, decreasing = TRUE)
```

> **Interpretation**: The genes are sorted by **log2 fold change**, which allows GSEA to identify pathways enriched in **upregulated** or **downregulated** genes.

---

### **Code: Running GSEA Analysis**  
We perform **GSEA using MSigDB hallmark gene sets (H)**.

```r
# Run GSEA analysis
```

> **Interpretation**: The output lists significantly enriched **biological pathways**, including:  
- **NES (Normalized Enrichment Score)**: Strength of pathway enrichment.  
- **FDR q-value**: Adjusted p-value indicating statistical significance.  

---

### **Code: Dot Plot of Enriched Pathways**  
We visualize the top **significantly enriched pathways**.

```r
# Filter pathways with NES ≥ 2 or ≤ -2 and FDR q-value ≤ 0.05
significant_pathways <- gsea_results@result[gsea_results@result$NES >= 2 | gsea_results@result$NES <= -2, ]
significant_pathways <- significant_pathways[significant_pathways$p.adjust <= 0.05, ]

# Create dot plot
ggplot(significant_pathways, aes(x = NES, y = reorder(Description, NES), size = setSize, color = p.adjust)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(low = "red", high = "blue") +
  labs(title = "GSEA: Enriched Pathways in WT vs. KO", x = "Normalized Enrichment Score (NES)", y = "Pathways") +
  theme_minimal()

# Save plot
ggsave("figures/GSEA_dotplot.pdf", width = 10, height = 8, dpi = 300)
```

> **Interpretation**:  
- **X-axis (NES)**: Shows pathways enriched in **upregulated** (positive NES) or **downregulated** (negative NES) genes.  
- **Dot size**: Number of genes involved in each pathway.  
- **Color gradient**: FDR significance (red = high, blue = low).  

---
