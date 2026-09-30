# scRNAseq-GBM-reanalysis

## Requirements

### 00_LoadData.R

Loads the raw count matrix and creates the unfiltered Seurat object (`gbm`) used for downstream analysis.

Saves the raw Seurat object as `data/gbm_raw_seurat.rds` so that later steps can load it directly without rebuilding the object from the count matrix.

Tested with:

- R 4.5.1
- Seurat 5.5.1
- data.table 1.18.4

### 01_QC_DoubletDetection.R

Calculates and visualizes QC metrics, confirms the existing QC filtering, detects potential doublets using scDblFinder, and removes predicted doublets.

Saves the QC-passed singlet object as `data/gbm_QC_singlets.rds` so that downstream analysis can start directly from the completed QC step.

Tested with:

- R 4.5.1
- Seurat 5.5.1
- ggplot2 4.0.3
- SingleCellExperiment 1.32.0
- scDblFinder 1.24.10

### 02_Normalization_HVG_Scaling.R

Loads the QC-passed singlet Seurat object, performs LogNormalize normalization, identifies 2,000 highly variable features, visualizes the top variable genes, and scales the expression data.

Saves the processed Seurat object as `data/gbm_normalized_HVG_scaled.rds` for use in the next downstream analysis step.

Tested with:

- R 4.5.1
- Seurat 5.5.1 

### 03_Dimensionality Reduction_Clustering_and_UMAP.R

Loads the normalized and scaled Seurat object, runs PCA, and generates a UMAP of the uncorrected data to visualise sample-driven structure (`group.by = "orig.ident"`). Harmony is then used to integrate the samples, the neighborhood graph is rebuilt on the Harmony embedding (`dims = 1:12`), and clusters are identified with `FindClusters(resolution = 0.5)`. A post-Harmony UMAP is produced so the before/after integration can be compared directly.

Saves the integrated, clustered Seurat object as `data/gbm_harmony_clustered.rds`.

Tested with:

- R 4.5.1
- Seurat 5.5.1
- harmony 2.0.5

### 04_CellType_Annotation.R

Loads the Harmony-clustered Seurat object and identifies marker genes for every cluster with `FindAllMarkers(only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)`. The top 10 markers per cluster are selected by `avg_log2FC` and exported, then checked against a canonical lineage marker panel (myeloid, T cell, oligodendrocyte, endothelial, fibroblast, pericyte, glial, astrocytic and cycling programs) using `DotPlot()`.

Broad cell types are assigned by combining the cluster-specific DE markers with the canonical marker evidence, and clusters without convincing normal-lineage markers are conservatively retained as `Candidate neoplastic`. Saves the annotated object as `data/gbm_broad_celltype_annotated.rds`.

`presto` is used as an optional accelerator for the Wilcoxon test behind `FindMarkers()`/`FindAllMarkers()`. It is not required: if it is not installed, comment out `library(presto)` and Seurat falls back to its own Wilcoxon implementation.

Tested with:

- R 4.5.1
- Seurat 5.5.1
- presto 1.1.0
- ggplot2 4.0.3

### 05_Neoplastic_Cell_Analysis.R

Subsets the candidate neoplastic cells from the broad-cell-type annotated object, recalculates variable features and PCA within the neoplastic compartment, selects informative PCs, and runs Harmony integration across patients. The neighborhood graph, neoplastic subclusters (`resolution = 0.5`, `dims = 1:10`), and a post-Harmony UMAP are then generated.

Cluster markers are identified with `FindAllMarkers()` and the subclusters are assigned manual Neftel states (AC / MES / NPC). Differential expression between `Recurrence` and `Initial` is computed separately within each state (`FindMarkers()`, `min.pct = 0.10`, `logfc.threshold = 0.25`), and the resulting up/down gene sets are tested for Biological Process enrichment with `enrichGO()` (`org.Hs.eg.db`, BH adjustment, `pvalueCutoff = 0.05`, `qvalueCutoff = 0.05`).

Saves the annotated Seurat objects as `data/gbm_neoplastic_clustered.rds`, `data/gbm_neoplastic_annotated.rds`, and `data/gbm_neoplastic_final.rds`, together with the marker, differential expression, and GO enrichment tables in `data/`.

Tested with:

- R 4.6.1
- Seurat 5.5.1
- harmony 2.0.5
- presto 1.1.0
- clusterProfiler 4.20.0
- org.Hs.eg.db 3.23.1
- ggplot2 4.0.3
- patchwork 1.3.2
- dplyr 1.2.1