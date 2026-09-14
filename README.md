# scRNAseq-GBM-reanalysis

## Requirements

### 00_LoadData.R

Loads the raw count matrix and creates the unfiltered Seurat object (`gbm`) used for downstream analysis.

Tested with:
- R 4.5.1
- Seurat 5.5.1
- data.table 1.18.4

### 01_QC_DoubletDetection.R

Calculates and visualizes QC metrics, confirms the existing QC filtering, detects potential doublets using scDblFinder, and removes predicted doublets.

Tested with:
- R 4.5.1
- Seurat 5.5.1
- ggplot2 4.0.3
- SingleCellExperiment 1.32.0
- scDblFinder 1.24.10