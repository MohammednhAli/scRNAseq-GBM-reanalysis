# ============================================================
# GBM scRNA-seq ANALYSIS
# Dataset: GSE229779
#
# Script: 02_Normalization_HVG_Scaling.R
#
# Purpose:
#   1. Load the QC-passed singlet Seurat object
#   2. Normalize gene-expression counts
#   3. Identify highly variable features (HVGs)
#   4. Visualize the most variable genes
#   5. Scale gene-expression values
# ============================================================


# ============================================================
# 1. INSTALL PACKAGE IF NEEDED
# ============================================================

# Uncomment this line only if Seurat is not installed.

# install.packages("Seurat")


# ============================================================
# 2. LOAD PACKAGE
# ============================================================

# Seurat provides the functions used for normalization,
# variable-feature selection, visualization and scaling.

library(Seurat)


# ============================================================
# 3. LOAD QC-PASSED SEURAT OBJECT
# ============================================================

# Load the final QC-passed Seurat object produced by
# 01_QC_DoubletDetection.R.
#
# This object contains the retained singlet cells after
# QC assessment and doublet removal.
#
# readRDS() loads the saved object directly, avoiding the need
# to repeat data loading, QC and doublet detection.

gbm <- readRDS(
  "data/gbm_QC_singlets.rds"
)


# Check the loaded Seurat object.

gbm

dim(gbm)


# ============================================================
# 4. NORMALIZATION
# ============================================================

# NormalizeData() normalizes gene-expression counts between cells.
#
# The "LogNormalize" method:
#
#   1. Divides each gene count by the total counts in that cell
#   2. Multiplies the result by a scale factor of 10,000
#   3. Log-transforms the normalized values
#
# This helps account for differences in sequencing depth
# between individual cells.

gbm <- NormalizeData(
  gbm,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)


# ============================================================
# 5. HIGHLY VARIABLE FEATURES (HVGs)
# ============================================================

# FindVariableFeatures() identifies genes showing high
# cell-to-cell variation.
#
# These genes are particularly informative for distinguishing
# biological differences between cells and are used in later
# dimensionality-reduction and clustering steps.
#
# The "vst" method is used to select the 2,000 most highly
# variable genes.

gbm <- FindVariableFeatures(
  gbm,
  selection.method = "vst",
  nfeatures = 2000
)


# ============================================================
# 6. VISUALIZE HIGHLY VARIABLE FEATURES
# ============================================================

# VariableFeatures() returns the genes identified as highly
# variable.
#
# Select the top 10 most variable genes for labeling.

top10 <- head(
  VariableFeatures(gbm),
  10
)


# VariableFeaturePlot() displays gene variability and highlights
# the genes selected as highly variable.

plot1 <- VariableFeaturePlot(gbm)


# LabelPoints() adds labels for the 10 most highly variable
# genes to the plot.
#
# repel = TRUE helps prevent overlapping gene labels.

LabelPoints(
  plot = plot1,
  points = top10,
  repel = TRUE
)


# ============================================================
# 7. SCALING
# ============================================================

# ScaleData() standardizes gene-expression values so that genes
# are on a comparable scale.
#
# For each gene, expression values are centered and scaled
# across cells.
#
# This prevents genes with larger expression values from
# dominating downstream analyses.
#
# By default, ScaleData() scales the variable features identified
# in the previous step.

gbm <- ScaleData(gbm)


# ============================================================
# 8. FINAL PROCESSED OBJECT
# ============================================================

# gbm now contains:
#
#   - QC-passed singlet cells
#   - normalized expression values
#   - 2,000 highly variable features
#   - scaled expression values
#
# The object is now ready for the next stage of the analysis.

gbm


# ============================================================
# 9. SAVE PROCESSED SEURAT OBJECT
# ============================================================

# Save the normalized, HVG-selected and scaled Seurat object
# as an RDS file.
#
# This allows the next analysis script to load the completed
# object directly without repeating these processing steps.

saveRDS(
  gbm,
  file = "data/gbm_normalized_HVG_scaled.rds"
)
