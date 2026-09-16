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