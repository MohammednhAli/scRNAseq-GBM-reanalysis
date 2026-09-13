# scRNAseq-GBM-reanalysis

## Requirements

### 00_LoadData.R

Loads the raw count matrix and creates the unfiltered Seurat object (`gbm`) used for downstream analysis.

Tested with:
- R 4.5.1
- Seurat 5.5.1
- data.table 1.18.4